# library(testthat); library(scrapper); source("test-runAllNeighborSteps.R")

x <- t(as.matrix(iris[,1:4]))

test_that("runAllNeighborSteps gives the same results as the reference", {
    res <- runAllNeighborSteps(x, num.threads=2, collapse.search=FALSE)

    umap.ref <- runUmap(x)
    expect_identical(res$runUmap, umap.ref)
    umap.single <- runAllNeighborSteps(x, runTsne.args=NULL, clusterGraph.args=NULL)
    expect_identical(umap.single$runUmap, umap.ref)

    tsne.ref <- runTsne(x)
    expect_identical(res$runTsne, tsne.ref)
    tsne.single <- runAllNeighborSteps(x, runUmap.args=NULL, clusterGraph.args=NULL)
    expect_identical(tsne.single$runTsne, tsne.ref)

    graph <- buildSnnGraph(x)
    clustering <- clusterGraph(graph)
    expect_identical(res$clusterGraph, clustering)
    expect_null(res$buildSnnGraph)
    cluster.single <- runAllNeighborSteps(x, runUmap.args=NULL, runTsne.args=NULL)
    expect_identical(cluster.single$clusterGraph, clustering)
})

test_that("runAllNeighborSteps works with collapsed neighbors", {
    BNPARAM <- BiocNeighbors::KmknnParam()
    res <- runAllNeighborSteps(x, BNPARAM=BNPARAM, num.threads=2)
    res.collapse <- runAllNeighborSteps(x, collapse.search=TRUE, BNPARAM=BNPARAM, num.threads=2)
    expect_identical(res, res.collapse)
})

test_that("runAllNeighborSteps works when returning the graph", {
    BNPARAM <- BiocNeighbors::KmknnParam()
    out <- runAllNeighborSteps(x, BNPARAM=BNPARAM, return.graph=TRUE, num.threads=2)
    expect_identical(ncol(x), out$buildSnnGraph$vertices)
    expect_identical(ncol(x), length(out$clusterGraph$membership))
})

test_that("runAllNeighborSteps returns sensibly when nothing is requested", {
    res <- runAllNeighborSteps(x, runUmap.args=NULL, runTsne.args=NULL, clusterGraph.args=NULL)
    expect_type(res, "list")
    expect_identical(length(res), 0L)
})
