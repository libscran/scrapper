# library(testthat); library(scrapper); source("test-aggregateAcrossGenes.R")

library(Matrix)

set.seed(9999)

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

    avagg <- aggregateAcrossGenes(x, sets, average=TRUE)
    expect_identical(names(avagg), names(sets))
    for (s in names(sets)) {
        current <- sets[[s]]
        expect_equal(Matrix::colSums(x[current[[1]],,drop=FALSE] * current[[2]]) / sum(current[[2]]), avagg[[s]])
    }

    expect_error(aggregateAcrossGenes(x, list(list(1L, 2, 3))), "contain two vectors")
    expect_error(aggregateAcrossGenes(x, list(list(1L, numeric(0)))), "same length")
    expect_error(aggregateAcrossGenes(x, list(list(-1L, 2))), "out-of-range")

    # Duplicated indices are just ignored.
    sets2 <- lapply(sets,
        function(x) {
            first <- x[[1]]
            list(
                rep(first, 2),
                c(x[[2]], runif(length(first)))
            )
        }
    )
    expect_identical(aggregateAcrossGenes(x, sets2), agg)
})

test_that("aggregateAcrossGenes works for other vector types", {
    x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
    rownames(x) <- sprintf("GENE_%s", seq_len(nrow(x)))

    sets <- list(
       foo = rbinom(nrow(x), 1, 0.2) == 1,
       bar = sample(rownames(x), 20),
       stuff = list(sample(rownames(x), 50), runif(50))
    )
    agg <- aggregateAcrossGenes(x, sets)

    refsets <- sets
    refsets$foo <- which(sets$foo)
    refsets$bar <- match(sets$bar, rownames(x))
    refsets$stuff[[1]] <- match(sets$stuff[[1]], rownames(x))
    refagg <- aggregateAcrossGenes(x, refsets)

    expect_identical(agg, refagg)

    # Duplicated and missing names are just ignored, with or without weights.
    sets2 <- sets
    sets2$bar <- c("foo", sets$bar, "whee", sets$bar)
    sets2$stuff <- local({
        idx <- sets$stuff[[1]]
        list(
            c(rep(idx, 2), "absent"),
            c(sets$stuff[[2]], runif(length(idx)), pi)
        )
    })
    expect_identical(aggregateAcrossGenes(x, sets2), agg)
})
