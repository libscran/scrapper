# library(testthat); library(scrapper); source("test-clusterGraph.R")

data <- matrix(rnorm(10000), ncol=1000)
out <- buildSnnGraph(data)
ptr <- buildSnnGraph(data, as.pointer=TRUE)

test_that("clusterGraph works correctly for multilevel", {
    clust <- clusterGraph(out, method="multilevel")
    expect_true(all(clust$membership >= 1L | clust$membership <= ncol(data)))
    for (i in seq_along(clust$levels)) {
        lev <- clust$levels[[i]]
        expect_true(all(lev >= 1L | lev <= ncol(data)))
    }

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="multilevel")
    expect_identical(clust, clust2)

    # Works with the graph.
    g <- igraph::graph(out$edges, n = out$vertices)
    igraph::E(g)$weight <- out$weights
    clusti <- clusterGraph(g, method="multilevel")
    expect_identical(clust, clusti)
})

test_that("clusterGraph works correctly for Leiden", {
    clust <- clusterGraph(out, method="leiden")
    expect_true(all(clust$membership >= 1L | clust$membership <= ncol(data)))

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="leiden")
    expect_identical(clust, clust2)
})

test_that("clusterGraph works correctly for Walktrap", {
    clust <- clusterGraph(out, method="walktrap")
    expect_true(all(clust$membership >= 1L | clust$membership <= ncol(data)))

    # Something sensible happens with a pointer.
    clust2 <- clusterGraph(ptr, method="walktrap")
    expect_identical(clust, clust2)
})
