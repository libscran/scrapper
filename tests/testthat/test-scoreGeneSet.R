# library(testthat); library(scrapper); source("test-scoreGeneSet.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
normed <- normalizeCounts(x, size.factors=centerSizeFactors(Matrix::colSums(x)))

test_that("scoreGeneSet works more or less as expected", {
    pcs <- scoreGeneSet(normed, 1:50)
    expect_identical(length(pcs$scores), ncol(x))
    expect_identical(nrow(pcs$weights), 50L)

    # Now with blocking. 
    pcs <- scoreGeneSet(normed, block=sample(3, ncol(x), replace=TRUE), 1:50)
    expect_identical(length(pcs$scores), ncol(x))
    expect_identical(nrow(pcs$weights), 50L)
})
