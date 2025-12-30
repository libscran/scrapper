# library(testthat); library(scrapper); source("test-rna_quality_control.R")

set.seed(444555666)
library(Matrix)
x <- round(abs(rsparsematrix(5000, 1000, 0.1) * 100))

test_that("computeRnaQcMetrics works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))
    expect_identical(qc$subsets$Mito, Matrix::colSums(x[sub$Mito,]) / qc$sum)

    expect_error(computeRnaQcMetrics(x, list(-1)), "names")
    expect_error(computeRnaQcMetrics(x, list(foo=TRUE)), "number of rows")
    expect_error(computeRnaQcMetrics(x, list(foo=-1)), "out-of-range")
    expect_error(computeRnaQcMetrics(x, list(foo="foo")), "no row names")
})

test_that("suggestRnaQcThresholds works as expected", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    num.mads <- 1.5
    thresholds <- suggestRnaQcThresholds(qc, num.mads=num.mads)

    # Check the thresholds.
    expect_lt(thresholds$sum, median(qc$sum))
    expect_lt(thresholds$detected, median(qc$detected))

    sub <- qc$subsets[[1]]
    expect_gt(thresholds$subsets[[1]], median(sub))
    expect_equal(names(thresholds$subsets), "Mito")

    # Check the filter.
    expected <- qc$sum >= thresholds$sum & qc$detected >= thresholds$detected & qc$subsets[[1]] <= thresholds$subsets[[1]]
    observed <- filterRnaQcMetrics(thresholds, qc)
    expect_identical(expected, observed)

    expect_error(filterRnaQcMetrics(thresholds, qc, block=1:10), "should be set to NULL")
})

test_that("suggestRnaQcThresholds works as expected with blocking", { 
    sub <- list(Mito=rbinom(nrow(x), 1, 0.1) > 0)
    qc <- computeRnaQcMetrics(x, sub)
    block <- sample(3, ncol(x), replace=TRUE)
    num.mads <- 1.5
    thresholds <- suggestRnaQcThresholds(qc, block=block, num.mads=num.mads)

    expect_identical(thresholds$block.ids, 1:3)

    # Check the thresholds.
    for (b in 1:3) {
        keep <- block == b
        expect_lt(thresholds$sum[[b]], median(qc$sum[keep]))
        expect_lt(thresholds$detected[[b]], median(qc$detected[keep]))
        expect_gt(thresholds$subsets[[1]][[b]], median(qc$subsets[[1]][keep]))
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

    expect_error(filterRnaQcMetrics(thresholds, qc), "expected 'block='")
    expect_error(filterRnaQcMetrics(thresholds, qc, block=1:10), "not present in 'thresholds'")
})
