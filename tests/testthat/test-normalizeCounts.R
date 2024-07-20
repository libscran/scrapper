# library(testthat); library(scrapper); source("test-normalizeCounts.R")

set.seed(8888)
library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
sf <- centerSizeFactors(colSums(x))

test_that("normalizeCounts gives the same results in all modes", {
    y <- normalizeCounts(x, sf)
    z <- normalizeCounts(beachmat::initializeCpp(x), sf)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, log=FALSE)
    z <- normalizeCounts(beachmat::initializeCpp(x), sf, log=FALSE)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, log.base=10)
    z <- normalizeCounts(beachmat::initializeCpp(x), sf, log.base=10)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, pseudo.count=3)
    z <- normalizeCounts(beachmat::initializeCpp(x), sf, pseudo.count=3)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, pseudo.count=3, preserve.sparsity=TRUE)
    z <- normalizeCounts(beachmat::initializeCpp(x), sf, pseudo.count=3, preserve.sparsity=TRUE)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))
})
