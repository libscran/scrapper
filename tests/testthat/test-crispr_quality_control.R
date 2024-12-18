# library(testthat); library(scrapper); source("test-crispr_quality_control.R")

library(Matrix)
x <- round(abs(rsparsematrix(25, 100, 0.1) * 100))
z <- as.matrix(x)

test_that("computeCrisprQcMetrics works as expected", { 
    qc <- computeCrisprQcMetrics(x)
    expect_identical(qc$sum, Matrix::colSums(x))
    expect_identical(qc$detected, Matrix::colSums(x > 0))

    top <- max.col(t(z))
    expect_equal(qc$max.value, z[cbind(top, seq_along(top))])

    # Can't compare top and top2 directly in case of ties.
    top2 <- qc$max.index
    expect_equal(qc$max.value, z[cbind(top2, seq_along(top2))])
})

test_that("suggestCrisprQcFilters works as expected", { 
    qc <- computeCrisprQcMetrics(x)
    thresholds <- suggestCrisprQcThresholds(qc)

    expected <- qc$max.value >= thresholds$max.value
    expected[is.na(expected)] <- TRUE
    observed <- filterCrisprQcMetrics(thresholds, qc)
    expect_identical(expected, observed)
})

test_that("suggestCrisprQcFilters works as expected with blocking", { 
    qc <- computeCrisprQcMetrics(x)
    block <- sample(3, ncol(x), replace=TRUE)
    thresholds <- suggestCrisprQcThresholds(qc, block=block)

    expected <- qc$max.value >= unname(thresholds$max.value[block])
    expected[is.na(expected)] <- TRUE
    observed <- filterCrisprQcMetrics(thresholds, qc, block=block)
    expect_identical(expected, observed)

    # Same filtering with just the last block.
    last <- block == 3
    last_observed <- filterCrisprQcMetrics(
        thresholds,
        lapply(qc, function(x) x[last]),
        block=block[last]
    )
    expect_identical(observed[last], last_observed)
})
