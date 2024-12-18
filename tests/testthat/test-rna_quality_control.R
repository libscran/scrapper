# library(testthat); library(scrapper); source("test-rna_quality_control.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))

test_that("computeRnaQcMetrics works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))
    expect_identical(qc$subsets$Mito, Matrix::colSums(x[sub$Mito,]) / qc$sum)
})

test_that("suggestRnaQcThresholds works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    thresholds <- suggestRnaQcThresholds(qc)

    # Check the thresholds.
    lsum <- log(qc$sum)
    expect_equal(thresholds$sum, exp(median(lsum) - 3 * mad(lsum)))
    ldet <- log(qc$detected)
    expect_equal(thresholds$detected, exp(median(ldet) - 3 * mad(ldet)))

    sub <- qc$subsets[[1]]
    expect_equal(thresholds$subsets[[1]], median(sub) + 3 * mad(sub))
    expect_equal(names(thresholds$subsets), "Mito")

    # Check the filter.
    expected <- qc$sum >= thresholds$sum & qc$detected >= thresholds$detected & qc$subsets[[1]] <= thresholds$subsets[[1]]
    observed <- filterRnaQcMetrics(thresholds, qc)
    expect_identical(expected, observed)
})

test_that("suggestRnaQcThresholds works as expected with blocking", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    block <- sample(3, ncol(x), replace=TRUE)
    thresholds <- suggestRnaQcThresholds(qc, block=block)

    # Check the thresholds.
    for (b in 1:3) {
        keep <- block == b
        lsum <- log(qc$sum[keep])
        expect_equal(thresholds$sum[[b]], exp(median(lsum) - 3 * mad(lsum)))
        ldet <- log(qc$detected[keep])
        expect_equal(thresholds$detected[[b]], exp(median(ldet) - 3 * mad(ldet)))
        sub <- qc$subsets[[1]][keep]
        expect_equal(thresholds$subsets[[1]][[b]], median(sub) + 3 * mad(sub))
    }

    # Check the filter.
    expected <- qc$sum >= thresholds$sum[block] & qc$detected >= thresholds$detected[block] & qc$subsets[[1]] <= thresholds$subsets[[1]][block]
    observed <- filterRnaQcMetrics(thresholds, qc, block=block)
    expect_identical(unname(expected), observed)

    # Same filtering with just the last block.
    last <- block == 3
    last_observed <- filterRnaQcMetrics(
        thresholds,
        list(sum=qc$sum[last], detected=qc$detected[last], subsets=list(qc$subsets[[1]][last])),
        block=block[last]
    )
    expect_identical(observed[last], last_observed)
})
