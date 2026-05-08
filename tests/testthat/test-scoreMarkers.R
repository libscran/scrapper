# library(testthat); library(scrapper); source("test-scoreMarkers.R")

set.seed(999)
library(Matrix)
x <- abs(rsparsematrix(1000, 100, 0.1))
x@x <- jitter(x@x)

check_summaries <- function(summary) {
    for (x in names(summary)) {
        current <- summary[[x]]
        expect_true(all(current$min <= current$max))
        expect_true(all(current$min <= current$mean + 1e-8)) # add some tolerance for numerical imprecision when averaging identical effects.
        expect_true(all(current$min <= current$median))
        expect_true(all(current$mean <= current$max + 1e-8))
        expect_true(all(current$median <= current$max))
        expect_true(all(current$min_rank <= 1000))
        expect_true(all(current$min_rank > 0))
    }
}

test_that("scoreMarkers works as expected for simple cases", {
    g <- sample(letters[1:4], ncol(x), replace=TRUE)
    out <- scoreMarkers(x, g)

    expect_equal(out$mean[,"a"], Matrix::rowMeans(x[,g=="a"]))
    expect_equal(out$detected[,"d"], Matrix::rowMeans(x[,g=="d"] > 0))

    check_summaries(out$cohens.d)
    check_summaries(out$auc)
    check_summaries(out$delta.mean)
    check_summaries(out$delta.detected)

    for (n in names(out$auc)) {
        collected <- unlist(out$auc[[n]][c("min", "mean", "median", "max")])
        expect_true(all(collected >= 0))
        expect_true(all(collected <= 1))
    }

    # Same results for multiple threads.
    pout <- scoreMarkers(x, g, num.threads=2)
    expect_identical(out, pout)

    # Works without the AUC and the groupwise stats.
    aout <- scoreMarkers(x, g, compute.auc=FALSE, compute.group.mean=FALSE, compute.group.detected=FALSE)
    expect_false("auc" %in% names(aout))
    expect_false("mean" %in% names(aout))
    expect_false("detected" %in% names(aout))
    expect_equal(out$delta.detected, aout$delta.detected)

    # Works without anything.
    empty <- scoreMarkers(x, g, compute.auc=FALSE, compute.cohens.d=FALSE, compute.delta.detected=FALSE, compute.delta.mean=FALSE)
    expect_identical(names(empty), c("mean", "detected", "nrow", "group.ids"))

    expect_error(scoreMarkers(SummarizedExperiment::SummarizedExperiment(x), g), "not supported")
    expect_error(scoreMarkers(x, NULL), "non-NULL")
    expect_error(scoreMarkers(x, 1), "length equal")
})

test_that("scoreMarkers works with blocking", {
    g <- sample(letters[1:4], ncol(x), replace=TRUE)
    b <- rep(1:3, length.out=ncol(x))
    out <- scoreMarkers(x, g, block=b, block.weight.policy="equal")

    bkeep <- g == "b"
    expect_equal(out$mean[,"b"], (Matrix::rowMeans(x[,bkeep & b == 1]) + Matrix::rowMeans(x[,bkeep & b == 2]) + Matrix::rowMeans(x[,bkeep & b == 3]))/3)
    ckeep <- g == "c"
    expect_equal(out$mean[,"c"], (Matrix::rowMeans(x[,ckeep & b == 1]) + Matrix::rowMeans(x[,ckeep & b == 2]) + Matrix::rowMeans(x[,ckeep & b == 3]))/3)

    check_summaries(out$cohens.d)
    check_summaries(out$auc)
    check_summaries(out$delta.mean)
    check_summaries(out$delta.detected)

    for (group in names(out$auc)) {
        collected <- unlist(out$auc[[group]][c("min", "mean", "median", "max")])
        expect_true(all(collected >= 0))
        expect_true(all(collected <= 1))
    }

    # Weighted mean, unweighted mean and quantile are the same when there are only two groups.
    g2 <- rep(LETTERS[1:2], length.out=ncol(x))
    b2 <- rep(1:2, each=ncol(x)/2)
    ref <- scoreMarkers(x, g2, block=b2)

    mout <- scoreMarkers(x, g2, block=b2, block.average.policy="mean", block.weight.policy="equal")
    expect_equal(ref$mean, mout$mean)
    expect_equal(ref$cohens.d$mean, mout$cohens.d$mean)
    expect_equal(ref$auc$median, mout$cohens.d$median)

    qout <- scoreMarkers(x, g2, block=b2, block.average.policy="quantile")
    expect_equal(ref$detected, qout$detected)
    expect_equal(ref$delta.mean$min, mout$delta.mean$min)
    expect_equal(ref$delta.detected$median, mout$delta.detected$median)
})

test_that("scoreMarkers works with quantile summaries", {
    g <- sample(4, ncol(x), replace=TRUE)
    res <- scoreMarkers(x, g, compute.summary.quantiles=c(0, 0.5, 1))
    for (i in seq_len(max(g))) {
        for (e in c("cohens.d", "auc", "delta.mean", "delta.detected")) {
            current <- res[[e]][[i]]
            expect_identical(current$min, current$quantile[["0"]])
            expect_identical(current$median, current$quantile[["0.5"]])
            expect_identical(current$max, current$quantile[["1"]])
        }
    }
})

test_that("scoreMarkers works with all summaries disabled", {
    g <- sample(4, ncol(x), replace=TRUE)
    res <- scoreMarkers(
        x, 
        g,
        compute.group.mean=FALSE,
        compute.group.detected=FALSE,
        compute.summary.min=FALSE,
        compute.summary.mean=FALSE,
        compute.summary.median=FALSE,
        compute.summary.max=FALSE,
        compute.summary.min.rank=FALSE
    )

    for (i in seq_len(max(g))) {
        for (e in c("cohens.d", "auc", "delta.mean", "delta.detected")) {
            expect_identical(ncol(res[[e]][[i]]), 0L)
            expect_identical(nrow(res[[e]][[i]]), nrow(x))
        }
    }
})

test_that("scoreMarkers works as expected for the full pairwise statistics", {
    g <- sample(4, ncol(x), replace=TRUE)
    full <- scoreMarkers(x, g, all.pairwise=TRUE)

    # Checking that we set the dimensions correctly.
    for (g1 in seq_len(4)) {
        for (g2 in seq_len(g1)) {
            if (g1 == g2) {
                expect_true(all(full$delta.mean[g1, g2, ] == 0))
                expect_true(all(full$auc[g1, g2, ] == 0))
            } else {
                expect_equal(full$delta.mean[g1, g2, ], -full$delta.mean[g2, g1, ])
                expect_equal(full$auc[g1, g2, ], 1 - full$auc[g2, g1, ])
            }
        }
    }

    expect_true(all(full$auc >= 0))
    expect_true(all(full$auc <= 1))

    # Works without AUCs and the groupwise stats.
    aout <- scoreMarkers(x, g, all.pairwise=TRUE, compute.auc=FALSE, compute.group.mean=FALSE, compute.group.detected=FALSE)
    expect_false("auc" %in% names(aout))
    expect_false("mean" %in% names(aout))
    expect_false("detected" %in% names(aout))
    expect_equal(full$cohens.d, aout$cohens.d)

    # Works without anything.
    empty <- scoreMarkers(x, g, compute.auc=FALSE, compute.cohens.d=FALSE, compute.delta.detected=FALSE, compute.delta.mean=FALSE)
    expect_identical(names(empty), c("mean", "detected", "nrow", "group.ids"))

    # Works with blocking.
    b <- rep(1:3, length.out=ncol(x))
    bout <- scoreMarkers(x, g, block=b, block.weight.policy="equal", all.pairwise=TRUE)
    sbout <- scoreMarkers(x, g, block=b, block.weight.policy="equal")
    expect_identical(bout$mean, sbout$mean)
})

test_that("scoreMarkers works as expected for the best pairwise statistics", {
    convert_to_best <- function(pairwise, n, bound) {
        all.labels <- dimnames(pairwise)[[1]]
        all.genes <- dimnames(pairwise)[[3]]
        output <- vector("list", length(all.labels))
        names(output) <- all.labels
        for (g1 in all.labels) {
            current <- vector("list", length(all.labels))
            names(current) <- all.labels
            for (g2 in all.labels) {
                if (g1 == g2) {
                    next
                }
                stats <- pairwise[g2, g1,] # remember, second dimension is the first group in the comparison.
                keep <- which(stats > bound)
                o <- order(stats[keep], decreasing=TRUE)
                indices <- keep[head(o, n)]
                current[[g2]] <- S4Vectors::DataFrame(index=indices, effect=stats[indices])
            }
            output[[g1]] <- current
        }
        output
    }

    g <- sample(4, ncol(x), replace=TRUE)
    full <- scoreMarkers(x, g, all.pairwise=TRUE)
    best <- scoreMarkers(x, g, all.pairwise=10)

    expect_equal(full$mean, best$mean)
    expect_equal(full$detected, best$detected)
    expect_identical(convert_to_best(full$cohens.d, 10, 0), best$cohens.d)
    expect_identical(convert_to_best(full$auc, 10, 0.5), best$auc)
    expect_identical(convert_to_best(full$delta.mean, 10, 0), best$delta.mean)
    expect_identical(convert_to_best(full$delta.detected, 10, 0), best$delta.detected)

    # Works without AUCs and the groupwise means.
    aout <- scoreMarkers(x, g, all.pairwise=10, compute.auc=FALSE, compute.group.mean=FALSE, compute.group.detected=FALSE)
    expect_false("auc" %in% names(aout))
    expect_false("mean" %in% names(aout))
    expect_false("detected" %in% names(aout))
    expect_equal(best$cohens.d, aout$cohens.d)

    # Same results if we only want the top indices.
    best.index <- scoreMarkers(x, g, all.pairwise=10, compute.auc=FALSE, compute.group=FALSE, compute.group.detected=FALSE, top.index.only=TRUE)
    expect_identical(names(best.index$cohens.d), names(best$cohens.d))
    for (n in names(best$cohens.d)) {
        expect_identical(names(best.index$cohens.d[[n]]), names(best$cohens.d[[n]]))
        for (n2 in names(best$cohens.d[[n]])) {
            expect_identical(best.index$cohens.d[[n]][[n2]], best$cohens.d[[n]][[n2]]$index)
        }
    }

    # Works with blocking.
    b <- rep(1:3, length.out=ncol(x))
    bout <- scoreMarkers(x, g, block=b, block.weight.policy="equal", all.pairwise=10)
    sbout <- scoreMarkers(x, g, block=b, block.weight.policy="equal")
    expect_identical(bout$mean, sbout$mean)
})

test_that("reportGroupMarkerStatistics works as expected", {
    g <- sample(letters[1:4], ncol(x), replace=TRUE)
    out <- scoreMarkers(x, g)

    df <- reportGroupMarkerStatistics(out, "b")
    expect_identical(length(df$mean), nrow(x))
    expect_identical(length(df$detected), nrow(x))
    expect_identical(length(df$cohens.d.min), nrow(x))
    expect_identical(length(df$auc.median), nrow(x))
    expect_identical(length(df$delta.mean.min.rank), nrow(x))
    expect_identical(length(df$delta.detected.max), nrow(x))

    df <- reportGroupMarkerStatistics(out, 3, include.mean=FALSE, include.detected=FALSE, effect.sizes=character(0), summaries=character(0))
    expect_identical(nrow(df), nrow(x))
    expect_identical(ncol(df), 0L)

    # Works with row names and missing statistics.
    rownames(x) <- sprintf("GENE_%i", seq_len(nrow(x)))
    out <- scoreMarkers(x, g, compute.auc=FALSE)
    df <- reportGroupMarkerStatistics(out, 1)
    expect_identical(rownames(df), rownames(x))
    expect_null(df$auc.median)
})

test_that("defaults work as expected", {
    def <- scoreMarkersDefaults()
    expect_true(all(names(def) %in% names(formals(scoreMarkers))))
    def <- scoreMarkersDefaults(all.pairwise=TRUE)
    expect_true(all(names(def) %in% names(formals(scoreMarkers))))
    def <- scoreMarkersDefaults(all.pairwise=100)
    expect_true(all(names(def) %in% names(formals(scoreMarkers))))
})
