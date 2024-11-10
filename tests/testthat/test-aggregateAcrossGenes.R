# library(testthat); library(scrapper); source("test-aggregateAcrossGenes.R")

library(Matrix)

test_that("aggregateAcrossGenes works for unweighted sets", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))

    sets <- list(
       foo = sample(nrow(x), 20),
       bar = sample(nrow(x), 10),
       whee = sample(nrow(x), 500)
    )

    agg <- aggregateAcrossGenes(x, sets)
    expect_identical(names(agg), names(sets))
    for (s in names(sets)) {
        expect_equal(Matrix::colSums(x[sets[[s]],,drop=FALSE]), agg[[s]])
    }

    agg <- aggregateAcrossGenes(x, sets, average=TRUE)
    expect_identical(names(agg), names(sets))
    for (s in names(sets)) {
        expect_equal(Matrix::colMeans(x[sets[[s]],,drop=FALSE]), agg[[s]])
    }

    expect_error(aggregateAcrossGenes(x, list(-1L)), "out-of-range")
})

test_that("aggregateAcrossGenes works for weighted sets", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))

    sets <- list(
       foo = list(sample(nrow(x), 20), runif(20)),
       bar = list(sample(nrow(x), 10), runif(10)),
       whee = list(sample(nrow(x), 500), runif(500))
    )

    agg <- aggregateAcrossGenes(x, sets)
    expect_identical(names(agg), names(sets))
    for (s in names(sets)) {
        current <- sets[[s]]
        expect_equal(Matrix::colSums(x[current[[1]],,drop=FALSE] * current[[2]]), agg[[s]])
    }

    agg <- aggregateAcrossGenes(x, sets, average=TRUE)
    expect_identical(names(agg), names(sets))
    for (s in names(sets)) {
        current <- sets[[s]]
        expect_equal(Matrix::colSums(x[current[[1]],,drop=FALSE] * current[[2]]) / sum(current[[2]]), agg[[s]])
    }

    expect_error(aggregateAcrossGenes(x, list(list(1L, 2, 3))), "length 2")
    expect_error(aggregateAcrossGenes(x, list(list(1L, numeric(0)))), "equal length")
})
