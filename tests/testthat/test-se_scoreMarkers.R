# library(testthat); library(scrapper); source("test-se_scoreMarkers.R")

library(S4Vectors)

test_that("finding an ordering statistic works as expected", {
    expect_identical(scrapper:::.findOrderBy(DataFrame(cohens.d.mean=1, auc.max=2), TRUE), "cohens.d.mean")
    expect_identical(scrapper:::.findOrderBy(DataFrame(cohens.d.max=1, auc.median=2), TRUE), "auc.median")
    expect_null(scrapper:::.findOrderBy(DataFrame(cohens.d.max=1, auc.median=2), FALSE))
    expect_null(scrapper:::.findOrderBy(DataFrame(cohens.d.max=1, auc.median=2), NULL))
    expect_null(scrapper:::.findOrderBy(DataFrame(stuff=2), TRUE))
    expect_identical(scrapper:::.findOrderBy(DataFrame(cohens.d.max=1, auc.median=2), "foobar"), "foobar")
})

set.seed(6900)
library(SummarizedExperiment)
mat <- matrix(rpois(1000, 10), ncol=20)
rownames(mat) <- sprintf("gene%s", seq_len(nrow(mat)))
se <- SummarizedExperiment(list(counts=mat))
se <- normalizeRnaCounts.se(se)
groups <- rep(LETTERS[1:4], 5)

test_that("scoreMarkers.se works as expected", {
    out <- scoreMarkers.se(se, groups)
    expect_identical(names(out), LETTERS[1:4])

    for (g in LETTERS[1:4]) {
        df <- out[[g]]
        expect_identical(nrow(df), nrow(se))
        expect_true(all(rownames(mat) %in% rownames(df)))

        expect_type(df$cohens.d.mean, "double")
        expect_type(df$auc.median, "double")
        expect_type(df$delta.mean.min.rank, "integer")
        expect_type(df$mean, "double")
        expect_type(df$detected, "double")

        expect_false(is.unsorted(-df$cohens.d.mean)) # i.e., the default order-by choice.
    }

    # Also works with non-character groupings.
    out2 <- scoreMarkers.se(se, as.integer(factor(groups)))
    expect_identical(names(out2), as.character(1:4))
    for (g in 1:4) {
        expect_identical(out2[[as.character(g)]], out[[LETTERS[g]]])
    }
})

test_that("scoreMarkers.se works with extra columns", {
    rowData(se)$symbol <- sprintf("SYMBOL-%s", seq_len(nrow(mat)))
    out <- scoreMarkers.se(se, groups, extra.columns="symbol")
    for (g in LETTERS[1:4]) {
        df <- out[[g]]
        expect_identical(df$symbol, rowData(se)$symbol[match(rownames(df), rownames(se))])
    }
})

test_that("scoreMarkers.se works with quantiles", {
    out <- scoreMarkers.se(se, groups, more.marker.args=list(compute.summary.quantiles=c(0, 0.5, 1)))
    for (g in LETTERS[1:4]) {
        df <- out[[g]]
        expect_equal(df$auc.quantile.0.5, df$auc.median)
        expect_equal(df$cohens.d.quantile.0, df$cohens.d.min)
        expect_equal(df$delta.detected.quantile.1, df$delta.detected.max)
    }
})

test_that("scoreMarkers.se behaves without any metrics", {
    out <- scoreMarkers.se(se, groups, more.marker.args=list(compute.group.mean=FALSE, compute.group.detected=FALSE, compute.cohens.d=FALSE))
    for (g in LETTERS[1:4]) {
        df <- out[[g]]
        expect_identical(nrow(df), nrow(se))
        expect_true(all(rownames(mat) %in% rownames(df)))

        expect_null(df$cohens.d.mean)
        expect_null(df$mean)
        expect_null(df$detected)
        expect_type(df$auc.mean, "double")

        expect_false(is.unsorted(-df$auc.mean)) # i.e., the next default order-by choice.
    }
})

test_that("scoreMarkers.se sorts by min.rank correctly", {
    out <- scoreMarkers.se(se, groups, order.by="cohens.d.min.rank")
    for (g in LETTERS[1:4]) {
        df <- out[[g]]
        expect_false(is.unsorted(df$cohens.d.min.rank))
    }
})

test_that("scoreMarkers.se can still make DFs when no statistics are computed", {
    out <- scoreMarkers.se(se, groups, more.marker.args=list(
        compute.cohens.d=FALSE,
        compute.auc=FALSE,
        compute.group.mean=FALSE,
        compute.group.detected=FALSE,
        compute.delta.mean=FALSE,
        compute.delta.detected=FALSE
    ))

    for (g in LETTERS[1:4]) {
        expect_identical(nrow(out[[g]]), nrow(se))
        expect_identical(rownames(out[[g]]), rownames(se))
        expect_identical(ncol(out[[g]]), 0L)
    }
})

test_that("previewMarkers works as expected", {
    out <- scoreMarkers.se(se, groups)

    preview <- previewMarkers(out[[1]])
    expect_identical(colnames(preview), c("mean", "detected", "lfc"))
    expect_identical(nrow(preview), 10L)

    order_preview <- previewMarkers(out[[1]], order.by=TRUE)
    expect_identical(colnames(order_preview), c("mean", "detected", "lfc", "cohens.d.mean"))
    expect_identical(nrow(order_preview), 10L)

    noninc_preview <- previewMarkers(out[[1]], columns=NULL, order.by=TRUE, include.order.by=FALSE)
    expect_identical(colnames(noninc_preview), character(0))
    expect_identical(nrow(noninc_preview), 10L)
    expect_identical(rownames(order_preview), rownames(noninc_preview))

    eff_preview <- previewMarkers(out[[1]], order.by=TRUE, include.order.by="effect")
    expect_identical(eff_preview$effect, order_preview$cohens.d.mean)
    expect_identical(nrow(eff_preview), 10L)
    expect_identical(rownames(order_preview), rownames(eff_preview))

    preview <- previewMarkers(out[[1]], rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]]))

    preview <- previewMarkers(out[[1]], order.by="auc.median")
    expect_identical(nrow(preview), 10L)

    preview <- previewMarkers(out[[1]], order.by="auc.median", rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]])[order(out[[1]]$auc.median, decreasing=TRUE)])

    preview <- previewMarkers(out[[1]], order.by="auc.min.rank", rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]])[order(out[[1]]$auc.min.rank)])
})
