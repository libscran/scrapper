# library(testthat); library(scrapper); source("test-combineFactors.R")

test_that("combineFactors works for multiple factors", {
    set.seed(10000)
    x <- sample(LETTERS[1:5], 100, replace=TRUE)
    y <- sample(3, 100, replace=TRUE)
    z <- sample(letters[1:4], 100, replace=TRUE)

    out <- combineFactors(list(x, y, z))
    expect_identical(out$levels[[1]][out$index], x)
    expect_identical(out$levels[[2]][out$index], y)
    expect_identical(out$levels[[3]][out$index], z)
})

test_that("combineFactors works for multiple factors with unused combinations", {
    set.seed(10002)
    x <- factor(sample(LETTERS[1:5], 10, replace=TRUE), LETTERS[1:5])
    y <- factor(sample(3, 10, replace=TRUE), 1:3)
    z <- factor(sample(letters[1:4], 10, replace=TRUE), letters[1:4])

    # Sanity check.
    out <- combineFactors(list(x, y, z))
    expect_lt(length(out$levels[[1]]), 60L)
    expect_identical(out$levels[[1]][out$index], as.character(x))
    expect_identical(out$levels[[2]][out$index], as.character(y))
    expect_identical(out$levels[[3]][out$index], as.character(z))

    out <- combineFactors(list(x, y, z), keep.unused=TRUE)
    expect_equal(length(out$levels[[1]]), 60L)
    expect_identical(out$levels[[1]][out$index], as.character(x))
    expect_identical(out$levels[[2]][out$index], as.character(y))
    expect_identical(out$levels[[3]][out$index], as.character(z))
})

test_that("combineFactors works for a single factor", {
    set.seed(10001)
    x <- sample(LETTERS[1:5], 100, replace=TRUE)
    out <- combineFactors(list(x))
    expect_identical(out$levels[[1]][out$index], x)

    # Drops unused elements.
    x <- factor(sample(LETTERS, 10, replace=TRUE), LETTERS)
    out <- combineFactors(list(x))
    expect_identical(out$levels[[1]][out$index], as.character(x))
    expect_lt(length(out$levels[[1]]), length(LETTERS))

    out <- combineFactors(list(x), keep.unused=TRUE)
    expect_identical(out$levels[[1]][out$index], as.character(x))
    expect_identical(out$levels[[1]], LETTERS)
})
