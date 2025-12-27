# library(testthat); library(scrapper); source("test-se_scoreMarkers.R")

test_that("guessing the dimnames works as expected", {
    mat <- matrix(0, 200, 10)
    rownames(mat) <- sprintf("foo%s", seq_len(nrow(mat)))
    colnames(mat) <- LETTERS[1:10]

    # Works on a matrix.
    out <- scrapper:::.guessDimnames(list(stuff=mat))
    expect_identical(out$nrow, 200L)
    expect_identical(out$rownames, rownames(mat))
    expect_identical(out$groups, LETTERS[1:10])

    # Works for a list of effect sizes.
    effects <- lapply(1:10, function(x) data.frame(X=seq_len(200), row.names=rownames(mat)))
    names(effects) <- LETTERS[1:10]
    out2 <- scrapper:::.guessDimnames(list(cohens.d=effects))
    expect_identical(out, out2)
})

test_that("finding an ordering statistic works as expected", {
    expect_identical(scrapper:::.findOrderBy(data.frame(cohens.d.mean=1, auc.max=2), TRUE), "cohens.d.mean")
    expect_identical(scrapper:::.findOrderBy(data.frame(cohens.d.max=1, auc.median=2), TRUE), "auc.median")
    expect_null(scrapper:::.findOrderBy(data.frame(cohens.d.max=1, auc.median=2), FALSE))
    expect_null(scrapper:::.findOrderBy(data.frame(cohens.d.max=1, auc.median=2), NULL))
    expect_null(scrapper:::.findOrderBy(data.frame(stuff=2), TRUE))
    expect_identical(scrapper:::.findOrderBy(data.frame(cohens.d.max=1, auc.median=2), "foobar"), "foobar")
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

test_that("previewMarkers works as expected", {
    out <- scoreMarkers.se(se, groups)

    preview <- previewMarkers(out[[1]], NULL)
    expect_identical(colnames(preview), c("mean", "detected", "lfc"))
    expect_identical(nrow(preview), 10L)

    preview <- previewMarkers(out[[1]], NULL, rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]]))

    preview <- previewMarkers(out[[1]], order.by="auc.median")
    expect_identical(nrow(preview), 10L)

    preview <- previewMarkers(out[[1]], order.by="auc.median", rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]])[order(out[[1]]$auc.median, decreasing=TRUE)])

    preview <- previewMarkers(out[[1]], order.by="auc.min.rank", rows=NULL)
    expect_identical(rownames(preview), rownames(out[[1]])[order(out[[1]]$auc.min.rank)])
})
