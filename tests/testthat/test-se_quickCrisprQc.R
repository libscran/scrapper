# library(testthat); library(scrapper); source("test-se_quickCrisprQc.R")

library(SummarizedExperiment)
mat <- matrix(rpois(1000, 1), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("quickCrisprQc.se works as expected", {
    out <- quickCrisprQc.se(se)
    expect_type(out$sum, "double")
    expect_type(out$detected, "integer")
    expect_type(out$max.value, "double")
    expect_type(out$max.index, "integer")
    expect_type(out$keep, "logical")
    expect_type(metadata(out)$qc$thresholds, "list")

    out2 <- quickCrisprQc.se(se, output.prefix="WHEE.")
    expect_identical(out2$WHEE.sum, out$sum)
    expect_identical(out2$WHEE.detected, out$detected)
    expect_identical(out2$WHEE.max.value, out$max.value)
    expect_identical(out2$WHEE.max.value, out$max.value)
    expect_identical(out2$WHEE.keep, out$keep)

    out4 <- quickCrisprQc.se(se, meta.name=NULL)
    expect_null(metadata(out4)$qc)
})

test_that("quickCrisprQc.se works with filtering", {
    set.seed(239)
    mat <- matrix(rpois(10000, 1), ncol=100)
    se <- SummarizedExperiment(list(counts=mat))
    out <- quickCrisprQc.se(se, more.suggest.args=list(num.mads=0))
    expect_false(all(out$keep))

    filtered <- quickCrisprQc.se(se, filter.cells=TRUE, more.suggest.args=list(num.mads=0))
    expect_lt(ncol(filtered), ncol(se))
    expect_equal(ncol(filtered), sum(out$keep))
    expect_true(all(filtered$keep))
})

test_that("quickCrisprQc.se works with custom thresholds", {
    # Auto-fills the information that wasn't provided.
    fixed <- list()
    out <- quickCrisprQc.se(se, thresholds=fixed)
    expect_identical(length(metadata(out)$qc$thresholds$max.value), 1L)

    # Used verbatim (some hacking required for renamed subsets).
    thresholds <- list(max.value=pi)
    out <- quickCrisprQc.se(se, thresholds=I(thresholds))
    expect_identical(metadata(out)$qc$thresholds, thresholds)
})

test_that("quickCrisprQc.se overwrites existing entries", {
    copy <- se
    copy$max.value <- "A"
    copy$max.index <- "B" 

    copy <- quickCrisprQc.se(copy)
    expect_type(copy$max.value, "double")
    expect_type(copy$max.index, "integer")
})
