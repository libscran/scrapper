# library(testthat); library(scrapper); source("test-utils.R")

test_that("transformFactor behaves with typical arguments", {
    chosen <- sample(length(LETTERS), 50, replace=TRUE)
    uchosen <- sort(unique(chosen))
    idx <- match(chosen, uchosen) - 1L

    out <- scrapper:::.transformFactor(LETTERS[chosen])
    expect_identical(out$index, idx)
    expect_identical(out$names, LETTERS[uchosen])

    # Respects the input type.
    out <- scrapper:::.transformFactor(chosen * 2L)
    expect_identical(out$index, idx)
    expect_identical(out$names, uchosen * 2L)
    expect_type(out$names, "integer")

    out <- scrapper:::.transformFactor(chosen / 2)
    expect_identical(out$index, idx)
    expect_identical(out$names, uchosen / 2)
    expect_type(out$names, "double")
})

test_that("transformFactors behaves with factors", {
    chosen <- sample(5, 20, replace=TRUE)
    out <- scrapper:::.transformFactor(factor(LETTERS[chosen], rev(LETTERS)))
    expect_identical(out$index, 26L - chosen)
    expect_identical(out$names, factor(rev(LETTERS), rev(LETTERS)))

    # Even if they're ordered.
    chosen <- sample(13, 50, replace=TRUE) * 2L
    out <- scrapper:::.transformFactor(ordered(letters[chosen], letters))
    expect_identical(out$index, chosen - 1L)
    expect_identical(out$names, ordered(letters, letters))
})
