# library(testthat); library(scrapper); source("test-buildSnnGraph.R")

test_that("buildSnnGraph works correctly", {
    data <- matrix(rnorm(10000), ncol=1000)
    out <- buildSnnGraph(data)
    expect_identical(length(out$edges), 2L * length(out$weights))  
    expect_true(all(out$weights > 0))
    expect_true(all(out$edges >= 1 & out$edges <= ncol(data)))

    # Same results when given an index, and throwing in some parallelization.
    idx <- BiocNeighbors::buildIndex(data, transposed=TRUE, BNPARAM=BiocNeighbors::AnnoyParam())
    out2 <- buildSnnGraph(idx, num.threads=2)
    expect_identical(out, out2)

    # Same results when given a matrix of indices.
    res <- BiocNeighbors::findKNN(idx, k=10, get.distance=FALSE, get.index="transposed", BNPARAM=BiocNeighbors::AnnoyParam())
    out3 <- buildSnnGraph(res) 
    expect_identical(out, out3)

    # Something sensible happens with a pointer.
    out <- buildSnnGraph(data, as.pointer=TRUE)
    expect_type(out, "externalptr")
})
