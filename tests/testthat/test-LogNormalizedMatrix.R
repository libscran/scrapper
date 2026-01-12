# library(testthat); library(scrapper); source("test-LogNormalizedMatrix.R")

test_that("LogNormalizedMatrix works for dense matrices", {
    mat <- matrix(rpois(1000, lambda=2), ncol=10)
    sf <- centerSizeFactors(colSums(mat))

    norm <- LogNormalizedMatrix(mat, sf)
    expect_s4_class(norm, "LogNormalizedMatrix")
    expect_identical(dim(norm), dim(mat))
    expect_false(is_sparse(norm))
    expect_output(show(norm), "LogNormalizedMatrix")

    ref <- log1p(t(t(mat) / sf)) / log(2)
    expect_equal(as.matrix(norm), ref)
    expect_equal(extract_array(norm, list(1:10, 1:10)), ref[1:10,1:10])

    norm <- LogNormalizedMatrix(mat, sf, pseudo.count=2)
    ref <- log2(t(t(mat) / sf) + 2)
    expect_equal(as.matrix(norm), ref)

    norm <- LogNormalizedMatrix(mat, sf, log.base=10)
    ref <- log10(t(t(mat) / sf) + 1)
    expect_equal(as.matrix(norm), ref)

    dimnames(mat) <- list(paste0("GENE-", 1:100), LETTERS[1:10])
    norm <- LogNormalizedMatrix(mat, sf)
    expect_identical(dimnames(mat), dimnames(norm))

    expect_error(LogNormalizedMatrix(mat, 1), "equal number of columns")
})

library(Matrix)
test_that("LogNormalizedMatrix works for sparse matrices", {
    smat <- abs(rsparsematrix(50, 20, density=0.21))
    sf <- centerSizeFactors(colSums(smat))

    norm <- LogNormalizedMatrix(smat, sf)
    expect_true(is_sparse(norm))

    ref <- as.matrix(log1p(t(t(smat) / sf)) / log(2))
    expect_equal(extract_array(norm, list(1:10, 1:10)), ref[1:10,1:10])

    expect_equal({
        tmp <- as.matrix(norm)
        dimnames(tmp) <- NULL
        tmp
    }, as.matrix(ref))

    extracted <- extract_sparse_array(norm, list(NULL, NULL))
    expect_s4_class(extracted, "SVT_SparseMatrix")
    expect_equal(as.matrix(extracted), ref)

    extracted <- extract_sparse_array(norm, list(20:50, 10:20))
    expect_s4_class(extracted, "SVT_SparseMatrix")
    expect_equal(as.matrix(extracted), ref[20:50,10:20])
})

test_that("LogNormalizedMatrix works with initializeCpp", {
    mat <- matrix(rpois(1000, lambda=1), ncol=10)
    sf <- centerSizeFactors(colSums(mat))

    norm <- LogNormalizedMatrix(mat, sf)
    ref <- as.matrix(norm)
    ptr <- beachmat::initializeCpp(norm)
    expect_equal(ref, beachmat::tatami.realize(ptr, num.threads=1))

    # Also works with sparse matrices.
    smat <- as(mat, "dgCMatrix")
    snorm <- LogNormalizedMatrix(smat, sf)
    sptr <- beachmat::initializeCpp(snorm)
    expect_equal(as(snorm, "dgCMatrix"), beachmat::tatami.realize(sptr, num.threads=1))
})
