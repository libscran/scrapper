# library(testthat); library(scrapper); source("test-subsampleByNeighbors.R")

set.seed(10000)
x <- matrix(runif(10000), nrow=10)

test_that("subsampleByNeighbors works as expected", {
    keep <- subsampleByNeighbors(x, num.neighbors=5, min.remaining=2)
    expect_lt(length(keep), ncol(x))
    expect_gt(min(keep), 0)
    expect_lte(max(keep), ncol(x))

    keep2 <- subsampleByNeighbors(x, num.neighbors=10, min.remaining=2)
    expect_lt(length(keep2), length(keep))
    expect_gt(min(keep2), 0)
    expect_lte(max(keep2), ncol(x))

    # Same results when given some nearest neighbor results.
    res <- BiocNeighbors::findKNN(x, k=10, transposed=TRUE, BNPARAM=BiocNeighbors::AnnoyParam(), get.index="transposed", get.distance="transposed")
    keep.alt <- subsampleByNeighbors(res, num.neighbors=10, min.remaining=2)
    expect_identical(keep2, keep.alt)

    expect_error(subsampleByNeighbors(x, num.neighbors=10, min.remaining=1000), "should not be greater")
})
