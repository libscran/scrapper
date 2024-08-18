# library(testthat); library(scrapper); source("test-runUmap.R")

library(BiocNeighbors)
x <- t(as.matrix(iris[,1:4]))

test_that("runUmap works in basic mode", {
    embed <- runUmap(x)
    expect_identical(nrow(embed), ncol(x))
    expect_identical(ncol(embed), 2L)

    again <- runUmap(x)
    expect_identical(embed, again) # check it's reproducible

    alt <- runUmap(x, num.neighbors=20)
    expect_identical(dim(alt), dim(embed))
    expect_false(identical(alt, embed)) # check that it has an effect

    res <- findKNN(x, transposed=TRUE, k=15, get.distance="transposed", get.index="transposed", BNPARAM=AnnoyParam())
    nnin <- runUmap(res)
    expect_identical(nnin, embed)
})
