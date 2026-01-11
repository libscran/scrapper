# library(testthat); library(scrapper); source("test-correctMnn.R")

set.seed(99991)
test_that("correctMnn works correctly in simple cases", {
    x <- matrix(rnorm(10000), nrow=10)
    b <- rep(1:3, c(500, 300, 200))
    x[,b==2] <- x[,b==2] + 3
    x[,b==3] <- x[,b==3] + 5

    means <- vapply(split(colMeans(x), b), mean, 0)
    expect_gt(means[2] - means[1], 2)
    expect_gt(means[3] - means[1], 4)

    corrected <- correctMnn(x, b)
    expect_identical(dim(corrected$corrected), dim(x))

    means <- vapply(split(colMeans(corrected$corrected), b), mean, 0)
    expect_lt(abs(means[2] - means[1]), 1) # means are within 1 SD of the truth.
    expect_lt(abs(means[3] - means[1]), 1)

    expect_error(correctMnn(SummarizedExperiment::SummarizedExperiment(x)), "not supported")
})
