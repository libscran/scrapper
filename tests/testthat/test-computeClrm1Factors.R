# library(testthat); library(scrapper); source("test-computeClrm1Factors.R")

set.seed(10000)
library(Matrix)
x <- abs(rsparsematrix(20, 1000, 0.1) * 10)

test_that("computeClrm1Factors works against a reference", {
    ref <- expm1(Matrix::colMeans(log1p(x[Matrix::rowSums(x)>0,,drop=FALSE])))
    obs <- computeClrm1Factors(x)
    expect_equal(ref, obs)

    par <- computeClrm1Factors(x, num.threads=2)
    expect_equal(obs, par)
})
