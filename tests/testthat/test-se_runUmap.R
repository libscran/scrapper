# library(testthat); library(scrapper); source("test-se_runUmap.R")

library(SingleCellExperiment)
set.seed(7777)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA") <- matrix(runif(1000), nrow=200)

test_that("runUmap.se works as expected", {
    out <- runUmap.se(sce)
    expect_true("UMAP" %in% reducedDimNames(out))
})
