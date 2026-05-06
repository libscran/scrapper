# library(testthat); library(scrapper); source("test-modelGeneVariances.R")

set.seed(121212)
library(Matrix)
x <- abs(rsparsematrix(1000, 100, 0.1) * 10)

test_that("modelGeneVariances works without blocking", {
    out <- modelGeneVariances(x)$statistics
    expect_equal(out$variances, apply(as.matrix(x), 1, var))
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

    expect_error(modelGeneVariances(SummarizedExperiment::SummarizedExperiment(x)), "not supported")

    nofit <- modelGeneVariances(x, fit.trend = FALSE)
    expect_identical(out[,c("means", "variances")], nofit$statistics)
})

test_that("modelGeneVariances works with blocking", {
    block <- sample(LETTERS[1:3], ncol(x), replace=TRUE)
    out <- modelGeneVariances(x, block, block.weight.policy="equal")

    expect_identical(out$block.ids, LETTERS[1:3])
    for (b in out$block.ids) {
        sub <- x[,block == b]
        current <- out$per.block[[b]]
        expect_equal(current$variances, apply(as.matrix(sub), 1, var))
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

    # Default is to use the mean.
    def <- modelGeneVariances(x, block)
    out2 <- modelGeneVariances(x, block, block.average.policy="mean")
    expect_identical(def, out2)

    # Quantile is the same as the mean with equal weights when there are two groups.
    block2 <- sample(LETTERS[4:5], ncol(x), replace=TRUE)
    def2 <- modelGeneVariances(x, block2, block.weight.policy="equal")
    out2 <- modelGeneVariances(x, block2, block.average.policy="quantile")
    expect_identical(def2, out2)

    # Works if we disable the average policy.
    noave <- modelGeneVariances(x, block, block.average.policy="none")
    expect_identical(dim(noave$statistics), c(nrow(x), 0L))
    expect_identical(dim(noave$statistics), c(nrow(x), 0L))
    expect_identical(out$per.block, noave$per.block)

    ref <- modelGeneVariances(x, block)
    nofit <- modelGeneVariances(x, block, fit.trend = FALSE)
    expect_identical(ref$block.ids, nofit$block.ids)
    expect_identical(ref$statistics[,c("means", "variances")], nofit$statistics)
    for (b in ref$block.ids) {
        expect_identical(ref$per.block[[b]][,c("means", "variances")], nofit$per.block[[b]])
    }

    # Works with other weight policies.
    noave <- modelGeneVariances(x, block, block.average.policy="none")
})

test_that("defaults work as expected", {
    def <- modelGeneVariancesDefaults()
    expect_true(all(names(def) %in% names(formals(modelGeneVariances))))
})
