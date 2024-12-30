# library(testthat); library(scrapper); source("test-runPca.R")

library(Matrix)
x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
normed <- normalizeCounts(x, size.factors=centerSizeFactors(Matrix::colSums(x)))

test_that("runPCA works more or less as expected", {
    pcs <- runPca(normed)

    rm <- rowMeans(pcs$components)
    expect_true(all(abs(rm) < 1e-8))
    expect_identical(length(rm), 25L)

    expect_equal(pcs$center, Matrix::rowMeans(normed))
    expect_null(pcs$scale)

    expect_false(is.unsorted(-pcs$variance.explained))
    rv <- apply(pcs$components, 1, var)
    expect_equal(rv, pcs$variance.explained)

    expect_identical(nrow(pcs$rotation), 1000L)
    expect_identical(ncol(pcs$rotation), 25L)

    # Works with scaling.
    spcs <- runPca(normed, scale=TRUE)
    expect_identical(length(spcs$scale), 1000L)
})

test_that("runPCA works with blocking", {
    block <- sample(3, ncol(x), replace=TRUE)
    pcs <- runPca(normed, block=block)

    rm <- rowMeans(pcs$components)
    expect_true(all(abs(rm) < 1e-8))
    expect_identical(length(rm), 25L)

    expect_equal(pcs$center[1,], Matrix::rowMeans(normed[,block==1]))
    expect_equal(pcs$center[2,], Matrix::rowMeans(normed[,block==2]))
    expect_equal(pcs$center[3,], Matrix::rowMeans(normed[,block==3]))

    # Variance isn't so easily computed this time, so we just check it's sorted.
    expect_false(is.unsorted(-pcs$variance.explained))

    expect_identical(nrow(pcs$rotation), 1000L)
    expect_identical(ncol(pcs$rotation), 25L)

    # Check that the PCs are different from the unblocked case.
    ref <- runPca(normed)
    expect_false(isTRUE(all.equal(ref$components, pcs$components)))
})

test_that("runPCA works with residual components", {
    block <- sample(3, ncol(x), replace=TRUE)
    pcs <- runPca(normed, block=block, components.from.residuals=TRUE)

    rm <- rowMeans(pcs$components)
    expect_true(all(abs(rm) < 1e-8))
    expect_identical(length(rm), 25L)

    expect_false(is.unsorted(-pcs$variance.explained))
    rv <- apply(pcs$components, 1, var)
    ratio <- rv / pcs$variance.explained
    expect_lt(diff(range(ratio)), 1e-8)

    expect_identical(nrow(pcs$rotation), 1000L)
    expect_identical(ncol(pcs$rotation), 25L)

    # Make sure it's different from the other options.
    ref1 <- runPca(normed)
    expect_false(isTRUE(all.equal(ref1$components, pcs$components)))
    ref2 <- runPca(normed, block=block)
    expect_false(isTRUE(all.equal(ref2$components, pcs$components)))
})
