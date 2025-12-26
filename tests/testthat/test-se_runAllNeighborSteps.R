# library(testthat); library(scrapper); source("test-se_runAllNeighborSteps.R")

library(SingleCellExperiment)
set.seed(6666)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA") <- matrix(runif(1000), nrow=200)

test_that("runAllNeighborSteps.se works as expected", {
    out <- runAllNeighborSteps.se(sce, num.threads=2)
    expect_true("TSNE" %in% reducedDimNames(out))
    expect_true("UMAP" %in% reducedDimNames(out))
    expect_true(is.factor(out$clusters))

    # Nullifying everything.
    null <- runAllNeighborSteps.se(sce, umap.output.name=NULL, tsne.output.name=NULL, cluster.output.name=NULL, num.threads=2)
    expect_false("TSNE" %in% reducedDimNames(null))
    expect_false("UMAP" %in% reducedDimNames(null))
    expect_null(null$clusters)

    # Checking that we can get the graph.
    wt.graph <- runAllNeighborSteps.se(sce, umap.output.name=NULL, tsne.output.name=NULL, build.graph.name="graph", num.threads=2)
    expect_false("TSNE" %in% reducedDimNames(wt.graph))
    expect_false("UMAP" %in% reducedDimNames(wt.graph))
    expect_true(is.factor(wt.graph$clusters))
    expect_type(metadata(wt.graph)$graph, "list")
})
