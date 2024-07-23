# library(testthat); library(scrapper); source("test-summarizeEffects.R")

set.seed(999)
library(Matrix)
x <- abs(rsparsematrix(1000, 100, 0.1))
x@x <- jitter(x@x)

test_that("summarizeEffects works as expected", {
    g <- sample(4, ncol(x), replace=TRUE)
    summ <- scoreMarkers(x, g)
    full <- scoreMarkers(x, g, all.pairwise=TRUE)
    expect_equal(summarizeEffects(full$cohens.d), summ$cohens.d)
    expect_equal(summarizeEffects(full$auc), summ$auc)
})
