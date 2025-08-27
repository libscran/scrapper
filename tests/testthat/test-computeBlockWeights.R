# library(testthat); library(scrapper); source("test-computeBlockWeights.R")

test_that("computeBlockWeights works as expected", {
    sizes <- c(1, 10, 100, 1000, 10000)
    expect_equal(computeBlockWeights(sizes), pmin(1, sizes/1000))
    expect_equal(computeBlockWeights(sizes, variable.block.weight=c(50, 500)), c(0, 0, 50/450, 1, 1))
    expect_equal(computeBlockWeights(sizes, block.weight.policy="equal"), rep(1, length(sizes)))
    expect_equal(computeBlockWeights(sizes, block.weight.policy="size"), sizes)
    expect_equal(computeBlockWeights(sizes, block.weight.policy="none"), sizes)
})
