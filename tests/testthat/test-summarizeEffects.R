# library(testthat); library(scrapper); source("test-summarizeEffects.R")

set.seed(999)
library(Matrix)
x <- abs(rsparsematrix(1000, 100, 0.1))
x@x <- jitter(x@x)

test_that("summarizeEffects works as expected", {
    g <- sample(4, ncol(x), replace=TRUE)
    summ <- scoreMarkers(x, g, min.rank.limit=nrow(x))
    full <- scoreMarkers(x, g, all.pairwise=TRUE)
    expect_equal(summarizeEffects(full$cohens.d), summ$cohens.d)
    expect_equal(summarizeEffects(full$auc), summ$auc)
})

test_that("summarizeEffects works with quantiles", {
    g <- sample(4, ncol(x), replace=TRUE)
    full <- scoreMarkers(x, g, all.pairwise=TRUE)
    res <- summarizeEffects(full$cohens.d, compute.summary.quantiles=c(0, 0.5, 1))

    for (i in seq_len(max(g))) {
        expect_identical(res[[i]]$median, res[[i]]$quantile[["0.5"]])
        expect_identical(res[[i]]$min, res[[i]]$quantile[["0"]])
        expect_identical(res[[i]]$max, res[[i]]$quantile[["1"]])
    }
})

test_that("summarizeEffects works with everything disabled", {
    g <- sample(4, ncol(x), replace=TRUE)
    full <- scoreMarkers(x, g, all.pairwise=TRUE)
    res <- summarizeEffects(
        full$cohens.d,
        compute.summary.min=FALSE,
        compute.summary.mean=FALSE,
        compute.summary.median=FALSE,
        compute.summary.max=FALSE,
        compute.summary.min.rank=FALSE
    )

    for (i in seq_len(max(g))) {
        expect_identical(ncol(res[[i]]), 0L)
        expect_identical(nrow(res[[i]]), nrow(x))
    }
})

