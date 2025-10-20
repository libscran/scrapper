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

test_that("sanitizeGeneSet works as expected", {
    expect_identical(scrapper:::.sanitizeGeneSet(5:50, n=100, names=NULL), 5:50)
    expect_identical(scrapper:::.sanitizeGeneSet(50:5, n=100, names=NULL), 5:50) # handles sorting
    expect_identical(scrapper:::.sanitizeGeneSet(c(5:50, 10:30), n=100, names=NULL), 5:50) # removes duplicates
    expect_error(scrapper:::.sanitizeGeneSet(5:50, n=10, names=NULL), "out-of-range")
    expect_error(scrapper:::.sanitizeGeneSet(NA_integer_, n=10, names=NULL), "out-of-range")

    x <- rbinom(100, 1, 0.5) == 1L
    expect_identical(scrapper:::.sanitizeGeneSet(x, n=100, names=NULL), which(x))
    expect_error(scrapper:::.sanitizeGeneSet(FALSE, n=10, names=NULL), "length of")

    expect_identical(scrapper:::.sanitizeGeneSet(LETTERS, n=100, names=LETTERS), 1:26)
    expect_identical(scrapper:::.sanitizeGeneSet(rev(LETTERS), n=100, names=LETTERS), 1:26) # handles sorting
    expect_identical(scrapper:::.sanitizeGeneSet(rep(LETTERS, each=2), n=100, names=LETTERS), 1:26)  # removes duplicates
    expect_error(scrapper:::.sanitizeGeneSet(LETTERS, n=10, names=NULL), "all elements")
})
