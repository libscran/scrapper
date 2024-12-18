# library(testthat); library(scrapper); source("test-adt_quality_control.R")

set.seed(112233)
library(Matrix)
x <- round(abs(rsparsematrix(200, 1000, 0.1) * 100))

test_that("computeAdtQcMetrics works as expected", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeAdtQcMetrics(x, sub)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))
    expect_identical(qc$subsets$IgG, Matrix::colSums(x[sub$IgG,]))
})

test_that("suggestAdtQcThresholds works as expected", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeAdtQcMetrics(x, sub)
    num.mads <- 1.5
    thresholds <- suggestAdtQcThresholds(qc, num.mads=num.mads)

    # Check the thresholds.
    expect_lt(thresholds$detected, median(qc$detected))
    expect_gt(thresholds$subsets[[1]], median(qc$subsets[[1]]))
    expect_equal(names(thresholds$subsets), "IgG")

    # Check the filter.
    expected <- qc$detected >= thresholds$detected & qc$subsets[[1]] <= thresholds$subsets[[1]]
    observed <- filterAdtQcMetrics(thresholds, qc)
    expect_identical(expected, observed) 
})

test_that("suggestAdtQcThresholds works as expected with blocking", { 
    sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeAdtQcMetrics(x, sub)
    block <- sample(3, ncol(x), replace=TRUE)
    num.mads <- 1.5
    thresholds <- suggestAdtQcThresholds(qc, block=block, num.mads=num.mads)

    # Check the thresholds.
    for (b in 1:3) {
        keep <- block == b
        expect_lt(thresholds$detected[b], median(qc$detected[keep]))
        expect_gt(thresholds$subsets[[1]][b], median(qc$subsets[[1]][keep]))
    }

    # Check the filter.
    expected <- qc$detected >= thresholds$detected[block] & qc$subsets[[1]] <= thresholds$subsets[[1]][block]
    observed <- filterAdtQcMetrics(thresholds, qc, block=block)
    expect_identical(unname(expected), observed)

    # Same filtering with just the last block.
    last <- block == 3
    last_observed <- filterAdtQcMetrics(
        thresholds,
        list(sum=qc$sum[last], detected=qc$detected[last], subsets=list(qc$subsets[[1]][last])),
        block=block[last]
    )
    expect_identical(observed[last], last_observed)
})
