# library(scrapper); library(testthat); source("test-testEnrichment.R")

test_that("testEnrichment works as expected", {
    sets <- list(
        first=LETTERS[1:10],
        second=LETTERS[1:5 * 2],
        third=LETTERS[10:20]
    )

    out <- testEnrichment(x=LETTERS[1:5], sets, universe=LETTERS)
    expected <- phyper(c(5, 2, 0) - 1L, lengths(sets), 26 - lengths(sets), 5L, lower.tail=FALSE)
    expect_identical(rownames(out), c("first", "second", "third"))
    expect_equal(out$p.value, expected)

    # Works with log-transformation.
    lout <- testEnrichment(x=LETTERS[1:5], sets, universe=LETTERS, log=TRUE)
    expected <- phyper(c(5, 2, 0) - 1L, lengths(sets), 26 - lengths(sets), 5L, log=TRUE, lower.tail=FALSE)
    expect_equal(lout$p.value, expected)

    # Works with parallelization.
    pout <- testEnrichment(x=LETTERS[1:5], sets, universe=LETTERS, num.threads=2)
    expect_equal(pout, out)
})

test_that("testEnrichment works with a truncated universe", {
    sets <- list(LETTERS[1:10])
    out <- testEnrichment(x=LETTERS[1:5], sets, universe=setdiff(LETTERS, "G"))
    expected <- phyper(5 - 1L, 9, 25 - 9, 5, lower.tail=FALSE)
    expect_equal(out$p.value, expected)

    sets <- list(LETTERS[1:10])
    out <- testEnrichment(x=LETTERS[1:5], sets, universe=setdiff(LETTERS, "A"))
    expected <- phyper(4 - 1L, 9, 25 - 9, 4, lower.tail=FALSE)
    expect_equal(out$p.value, expected)
})

test_that("testEnrichment works with a NULL universe", {
    sets <- list(
        first=LETTERS[1:10],
        second=LETTERS[1:5 * 2],
        third=LETTERS[10:20]
    )

    ref <- testEnrichment(x=LETTERS[1:5], sets, universe=LETTERS[1:20])
    out <- testEnrichment(x=LETTERS[1:5], sets, universe=NULL)
    expect_identical(ref, out)

    ref <- testEnrichment(x=LETTERS[c(1:4, 22:26)], sets, universe=LETTERS[-21])
    out <- testEnrichment(x=LETTERS[c(1:4, 22:26)], sets, universe=NULL)
    expect_identical(ref, out)
})

test_that("testEnrichment works with an integer universe", {
    sets <- list(LETTERS[1:10])
    out <- testEnrichment(x=LETTERS[1:5], sets, universe=26)
    expected <- phyper(5 - 1L, 10, 26 - 10, 5, lower.tail=FALSE)
    expect_equal(out$p.value, expected)
})
