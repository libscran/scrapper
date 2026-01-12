# library(testthat); library(scrapper); source("test-normalizeCounts.R")

set.seed(8888)
library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
sf <- centerSizeFactors(Matrix::colSums(x))

test_that("normalizeCounts gives the same results in all modes", {
    y <- normalizeCounts(x, sf)
    expect_s4_class(y, "LogNormalizedMatrix")
    z <- normalizeCounts(beachmat::initializeCpp(x), sf)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, delayed=FALSE)
    expect_s4_class(y, "dgCMatrix")
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, log=FALSE)
    expect_s4_class(y, "DelayedMatrix")
    z <- normalizeCounts(beachmat::initializeCpp(x), sf, log=FALSE)
    expect_equal(y[1,], beachmat:::tatami_row(z, 1))
    expect_equal(y[,1], beachmat:::tatami_column(z, 1))

    y <- normalizeCounts(x, sf, log=FALSE, delayed=FALSE)
    expect_s4_class(y, "dgCMatrix")
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

    expect_error(normalizeCounts(SummarizedExperiment::SummarizedExperiment(x)), "not supported")
})
