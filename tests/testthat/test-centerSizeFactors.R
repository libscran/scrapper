# library(testthat); library(scrapper); source("test-centerSizeFactors.R")

set.seed(99999)

test_that("centerSizeFactors works as expected", {
    f <- runif(100)
    cf <- centerSizeFactors(f)
    expect_equal(mean(cf), 1)
    expect_true(abs(diff(range(cf/f))) < 1e-8)

    # Ignores weirdos.
    f[1] <- 0
    f[2] <- NA
    f[3] <- Inf
    f[4] <- -1
    cf <- centerSizeFactors(f)
    tmp <- centerSizeFactors(tail(f, -4))
    expect_equal(tail(cf, -4), tmp)
})

test_that("centerSizeFactors works with blocking", {
    f <- runif(100)
    b <- sample(3, 100, replace=TRUE)
    scaling <- vapply(split(f, b), mean, 0)

    cf <- centerSizeFactors(f, b)
    expect_equal(f / min(scaling), cf)

    cf2 <- centerSizeFactors(f, b, mode="lowest")
    expect_equal(cf, cf2)

    cf <- centerSizeFactors(f, b, mode="per-block")
    expect_equal(f / unname(scaling[b]), cf)

    expect_error(centerSizeFactors(f, b, mode="WHEE"), "unknown")
})

test_that("defaults work correctly", {
    def <- centerSizeFactorsDefaults()
    expect_true(all(names(def) %in% names(formals(centerSizeFactors))))
})
