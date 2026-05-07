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

    expect_error(sanitizeSizeFactors(f, handle.zero="error"), "zero")
    expect_error(sanitizeSizeFactors(f, handle.nan="error"), "NaN")
    expect_error(sanitizeSizeFactors(f, handle.infinite="error"), "infinite")
    expect_error(sanitizeSizeFactors(f, handle.negative="error"), "negative")

    expect_identical(f, sanitizeSizeFactors(f, handle.zero="ignore", handle.negative="ignore", handle.infinite="ignore", handle.nan="ignore"))
})

test_that("defaults work correctly", {
    def <- sanitizeSizeFactorsDefaults()
    expect_true(all(names(def) %in% names(formals(sanitizeSizeFactors))))
})
