# library(testthat); library(scrapper); source("test-runTsne.R")

library(BiocNeighbors)
x <- t(as.matrix(iris[,1:4]))

test_that("runTSNE works in basic mode", {
    embed <- runTsne(x)
    expect_identical(nrow(embed), ncol(x))
    expect_identical(ncol(embed), 2L)

    again <- runTsne(x)
    expect_identical(embed, again)  # check that it's reproducible.

    alt <- runTsne(x, perplexity=20)
    expect_identical(dim(alt), dim(embed))
    expect_false(identical(alt, embed)) # check that perplexity has an effect.

    res <- findKNN(x, transposed=TRUE, k=tsnePerplexityToNeighbors(30), get.distance="transposed", get.index="transposed", BNPARAM=AnnoyParam())
    nnin <- runTsne(res)
    expect_identical(nnin, embed)
})
