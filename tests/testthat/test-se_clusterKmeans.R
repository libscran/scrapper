# library(testthat); library(scrapper); source("test-se_clusterKmeans.R")

library(SingleCellExperiment)
set.seed(8888)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA") <- matrix(runif(1000), nrow=200)

test_that("clusterKmeans.se works as expected", {
    out <- clusterKmeans.se(sce, k=2)
    expect_identical(nlevels(out$clusters), 2L)

    out <- clusterKmeans.se(sce, k=10, meta.name="kmeans")
    expect_identical(nlevels(out$clusters), 10L)
    expect_identical(ncol(metadata(out)$kmeans$centers), 10L)
})
