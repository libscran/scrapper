# library(testthat); library(scrapper); source("test-runUmap.R")

library(BiocNeighbors)
x <- t(as.matrix(iris[,1:4]))
ref <- runUmap(x)

test_that("runUmap works in basic mode", {
    expect_identical(nrow(ref), ncol(x))
    expect_identical(ncol(ref), 2L)

    again <- runUmap(x)
    expect_identical(ref, again) # check it's reproducible

    alt <- runUmap(x, num.neighbors=20)
    expect_identical(dim(alt), dim(ref))
    expect_false(identical(alt, ref)) # check that it has an effect

    res <- findKNN(x, transposed=TRUE, k=15, get.distance="transposed", get.index="transposed", BNPARAM=AnnoyParam())
    nnin <- runUmap(res)
    expect_identical(nnin, ref)
})

test_that("runUmap works with initial coordinates", {
    expect_error(runUmap(x, initialize.method="none"), "expected initial coordinates")
    expect_error(runUmap(x, initialize.method="none", initial.coordinates=cbind(1,2)), "nrow")
    expect_error(runUmap(x, initialize.method="none", initial.coordinates=cbind(seq_len(ncol(x)))), "ncol")

    stuff <- matrix(rnorm(ncol(x) * 2), ncol=2)
    init.out <- runUmap(x, initialize.method="none", initial.coordinates=stuff)
    expect_false(identical(init.out, ref))
})

test_that("runUmap works with pre-specified a/b", {
    custom <- runUmap(x, a=2, b=1)
    expect_false(identical(custom, ref))

    expect_error(runUmap(x, a=1:2, b=1), "should be a numeric scalar")
    expect_error(runUmap(x, a=2, b=1:2), "should be a numeric scalar")
})

test_that("runUmap works with an explicit number of epochs", {
    custom <- runUmap(x, num.epochs=600)
    expect_false(identical(custom, ref))
    expect_error(runUmap(x, num.epochs=1:1000), "should be an integer scalar")
})
