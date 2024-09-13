# library(testthat); library(scrapper); source("test-clusterKmeans.R")

x <- t(as.matrix(iris[,1:4]))

test_that("clusterKmeans works in basic mode", {
    clustering <- clusterKmeans(x, 10)
    expect_identical(length(clustering$clusters), ncol(x))
    expect_identical(length(unique(clustering$clusters)), 10L)
    expect_identical(ncol(clustering$centers), 10L)

    # Randomness should be fully controlled.
    again <- clusterKmeans(x, 10)
    expect_identical(clustering, again)
})

test_that("clusterKmeans works with different methods", {
    clustering <- clusterKmeans(x, 5, init.method="random")
    expect_identical(length(clustering$clusters), ncol(x))
    expect_identical(length(unique(clustering$clusters)), 5L)
    expect_identical(ncol(clustering$centers), 5L)

    clustering <- clusterKmeans(x, 5, init.method="kmeans++", refine.method="lloyd")
    expect_identical(length(clustering$clusters), ncol(x))
    expect_identical(length(unique(clustering$clusters)), 5L)
    expect_identical(ncol(clustering$centers), 5L)
})
