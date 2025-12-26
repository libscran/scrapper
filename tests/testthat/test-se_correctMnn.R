# library(testthat); library(scrapper); source("test-se_correctMnn.R")

library(SingleCellExperiment)
set.seed(8888)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA") <- matrix(runif(1000), nrow=200)

test_that("correctMnn.se works as expected", {
    block <- rep(LETTERS[1:2], each=100)
    out <- correctMnn.se(sce, block=block)
    expect_true(is.matrix(reducedDim(out, "MNN")))

    out <- correctMnn.se(sce, block=block, delayed.transpose=TRUE)
    expect_s4_class(reducedDim(out, "MNN"), "DelayedMatrix")
})
