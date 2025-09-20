# library(testthat); library(scrapper); source("test-chooseHighlyVariableGenes.R")

test_that("chooseHighlyVariableGenes works", {
    set.seed(20000)

    stats <- runif(10000)
    out <- chooseHighlyVariableGenes(stats, top=2000)
    expect_equal(length(out), 2000)
    expect_gt(min(stats[out]), max(stats[-out]))

    out <- chooseHighlyVariableGenes(stats, top=Inf)
    expect_identical(out, seq_along(stats))

    out <- chooseHighlyVariableGenes(stats, larger=FALSE, bound=NULL)
    expect_equal(length(out), 4000)
    expect_lt(max(stats[out]), min(stats[-out]))

    out <- chooseHighlyVariableGenes(stats, bound=0.9)
    expect_lt(length(out), 2000)
    expect_gt(min(stats[out]), 0.9)

    out <- chooseHighlyVariableGenes(numeric(10000), keep.ties=TRUE, bound=NULL)
    expect_equal(out, seq_len(10000))
})
