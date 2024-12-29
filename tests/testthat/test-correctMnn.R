# library(testthat); library(scrapper); source("test-correctMnn.R")

set.seed(99991)
test_that("correctMnn works correctly in simple cases", {
    x <- matrix(rnorm(10000), nrow=10)
    b <- rep(1:3, c(500, 300, 200))
    x[,b==2] <- x[,b==2] + 3
    x[,b==3] <- x[,b==3] + 5

    means <- vapply(split(colMeans(x), b), mean, 0)
    expect_gt(means[2] - means[1], 2)
    expect_gt(means[3] - means[2], 1)

    corrected <- correctMnn(x, b)
    expect_identical(dim(corrected$corrected), dim(x))
    expect_identical(length(corrected$num.pairs), 2L)
    expect_identical(corrected$merge.order, as.character(1:3)) # largest block always have the most RSS, and the most MNN pairs.

    means <- vapply(split(colMeans(corrected$corrected), b), mean, 0)
    expect_lt(abs(means[2] - means[1]), 0.5)
    expect_lt(abs(means[3] - means[2]), 0.5)
})

set.seed(99992)
test_that("correctMnn works with a specified order", {
    x <- t(as.matrix(iris[,1:4]))
    b <- sample(1:2, ncol(x), replace=TRUE)
    corrected <- correctMnn(x, b, order=c("2", "1"))

    expect_identical(dim(corrected$corrected), dim(x))
    expect_identical(length(corrected$num.pairs), 1L)
    expect_identical(corrected$merge.order, as.character(2:1))

    # Actually has an effect.
    corrected2 <- correctMnn(x, b, order=c("1", "2"))
    expect_false(identical(corrected$corrected, corrected2$corrected))

    # Errors out correctly.
    expect_error(correctMnn(x, b, order=c("1", "3")), "out-of-range")
    expect_error(correctMnn(x, b, order=c("1", "1")), "duplicate")
    expect_error(correctMnn(x, b, order=c("1")), "number of batches")
})
