# library(testthat); library(scrapper); source("test-subsampleByNeighbors.R")

set.seed(10000)
x <- matrix(runif(10000), nrow=10)

test_that("subsampleByNeighbors works as expected", {
    keep <- subsampleByNeighbors(x, num.neighbors=5, min.remaining=2)
    expect_true(length(keep) < ncol(x))
    expect_true(min(keep) > 0)
    expect_true(max(keep) <= ncol(x))

    keep2 <- subsampleByNeighbors(x, num.neighbors=10, min.remaining=2)
    expect_true(length(keep2) < length(keep))
    expect_true(min(keep2) > 0)
    expect_true(max(keep2) <= ncol(x))

    expect_error(subsampleByNeighbors(x, num.neighbors=10, min.remaining=1000), "should not be greater")
})
