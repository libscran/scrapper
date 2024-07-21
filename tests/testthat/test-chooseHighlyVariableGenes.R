# library(testthat); library(scrapper); source("test-chooseHighlyVariableGenes.R")

test_that("chooseHighlyVariableGenes works", {
    set.seed(20000)

    stats <- runif(10000)
    out <- chooseHighlyVariableGenes(stats, top=2000)
    expect_equal(length(out), 2000)
    expect_gt(min(stats[out]), max(stats[-out]))

    out <- chooseHighlyVariableGenes(stats, larger=FALSE)
    expect_equal(length(out), 4000)
    expect_lt(max(stats[out]), min(stats[-out]))

    out <- chooseHighlyVariableGenes(numeric(10000), keep.ties=TRUE)
    expect_equal(out, seq_len(10000))
})
