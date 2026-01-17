# library(testthat); library(scrapper); source("test-se_clusterGraph.R")

library(SingleCellExperiment)
set.seed(9999)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA") <- matrix(runif(1000), nrow=200)

test_that("clusterGraph.se works as expected", {
    sce <- clusterGraph.se(sce)
    expect_true(is.factor(sce$clusters))

    sce2 <- clusterGraph.se(sce, resolution=2)
    expect_false(identical(sce$clusters, sce2$clusters))

    scel <- clusterGraph.se(sce, method="leiden")
    expect_false(identical(sce$clusters, scel$clusters))

    sce <- clusterGraph.se(sce, graph.name="graph", meta.name="cluster")
    expect_type(metadata(sce)$cluster, "list")
    expect_null(metadata(sce)$cluster$membership)
    expect_type(metadata(sce)$graph, "list")
})
