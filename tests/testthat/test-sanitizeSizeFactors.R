# library(testthat); library(scrapper); source("test-sanitizeSizeFactors.R")

set.seed(99999)

test_that("sanitizeSizeFactors works as expected", {
    f <- 2^rnorm(100)
    f[1] <- 0
    f[2] <- NA
    f[3] <- Inf
    f[4] <- -1

    cf <- sanitizeSizeFactors(f)
    expect_true(all(is.finite(cf) & cf > 0))
    expect_equal(tail(cf, -4), tail(f, -4))
})
