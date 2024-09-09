# library(testthat); library(scrapper); source("test-modelGeneVariances.R")

set.seed(121212)
library(MatrixGenerics)
library(Matrix)
x <- abs(rsparsematrix(1000, 100, 0.1) * 10)

test_that("modelGeneVariances works without blocking", {
    out <- modelGeneVariances(x)$statistics
    expect_equal(out$variances, rowVars(x))
    expect_equal(out$means, Matrix::rowMeans(x))
    expect_equal(out$variances - out$fitted, out$residuals)

    fit <- fitVarianceTrend(out$means, out$variances)
    expect_identical(out$fitted, fit$fitted)
    expect_identical(out$residuals, fit$residuals)

    # Responds to trend-fitting options. 
    out2 <- modelGeneVariances(x, span=0.5)$statistics
    expect_identical(out$variances, out2$variances)
    expect_identical(out$means, out2$means)
    expect_false(identical(out$fitted, out2$fitted))

    fit2 <- fitVarianceTrend(out2$means, out2$variances, span=0.5)
    expect_equal(out2$fitted, fit2$fitted)
    expect_equal(out2$residuals, fit2$residuals)
})

test_that("modelGeneVariances works with blocking", {
    block <- sample(LETTERS[1:3], ncol(x), replace=TRUE)
    out <- modelGeneVariances(x, block, block.weight.policy="equal")

    for (b in LETTERS[1:3]) {
        sub <- x[,block == b]
        current <- out$per.block[[b]]
        expect_equal(current$variances, rowVars(sub))
        expect_equal(current$means, Matrix::rowMeans(sub))
        expect_equal(current$variances - current$fitted, current$residuals)
    }

    average.mean <- rowMeans(do.call(cbind, lapply(out$per.block, function(x) x$means)))
    expect_equal(average.mean, out$statistics$means)
    average.variance <- rowMeans(do.call(cbind, lapply(out$per.block, function(x) x$variances)))
    expect_equal(average.variance, out$statistics$variances)
    average.fitted <- rowMeans(do.call(cbind, lapply(out$per.block, function(x) x$fitted)))
    expect_equal(average.fitted, out$statistics$fitted)
    average.residuals <- rowMeans(do.call(cbind, lapply(out$per.block, function(x) x$residuals)))
    expect_equal(average.residuals, out$statistics$residuals)
})
