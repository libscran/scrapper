# library(testthat); library(scrapper); source("test-clusterGraph.R")

data <- matrix(rnorm(10000), ncol=1000)
out <- buildSnnGraph(data)
ptr <- buildSnnGraph(data, as.pointer=TRUE)

test_that("clusterGraph works correctly for multilevel", {
    clust <- clusterGraph(out, method="multilevel")
    expect_identical(length(clust$membership), ncol(data))
    expect_lte(nlevels(clust$membership), ncol(data))
    expect_gte(nlevels(clust$membership), 1L)
    expect_false(anyNA(clust$membership))

    for (cl in clust$levels) {
        expect_identical(length(cl), ncol(data))
        expect_lte(nlevels(cl), ncol(data))
        expect_gte(nlevels(cl), 1L)
        expect_false(anyNA(cl))
    }

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="multilevel")
    expect_identical(clust, clust2)

    # Works with the graph.
    g <- igraph::make_undirected_graph(out$edges, n = out$vertices)
    igraph::E(g)$weight <- out$weights
    clusti <- clusterGraph(g, method="multilevel")
    expect_identical(clust, clusti)
})

test_that("clusterGraph works correctly for Leiden", {
    clust <- clusterGraph(out, method="leiden")
    expect_identical(length(clust$membership), ncol(data))
    expect_lte(nlevels(clust$membership), ncol(data))
    expect_gte(nlevels(clust$membership), 1L)
    expect_false(anyNA(clust$membership))

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="leiden")
    expect_identical(clust, clust2)
})

test_that("clusterGraph works correctly for Walktrap", {
    clust <- clusterGraph(out, method="walktrap")
    expect_identical(length(clust$membership), ncol(data))
    expect_lte(nlevels(clust$membership), ncol(data))
    expect_gte(nlevels(clust$membership), 1L)
    expect_false(anyNA(clust$membership))

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="walktrap")
    expect_identical(clust, clust2)
})
