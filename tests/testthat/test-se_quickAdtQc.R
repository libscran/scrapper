# library(testthat); library(scrapper); source("test-se_quickAdtQc.R")

library(SummarizedExperiment)
mat <- matrix(rpois(1000, 1), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("quickAdtQc.se works as expected", {
    out <- quickAdtQc.se(se, subsets=list())
    expect_identical(ncol(out), ncol(se))
    expect_type(out$sum, "double")
    expect_type(out$detected, "integer")
    expect_type(out$keep, "logical")
    expect_type(metadata(out)$qc$thresholds, "list")

    out <- quickAdtQc.se(se, subsets=list(igg=1:5))
    expect_type(out$subset.sum.igg, "double")

    out2 <- quickAdtQc.se(se, subsets=list(igg=1:5), output.prefix="WHEE.")
    expect_identical(out2$WHEE.sum, out$sum)
    expect_identical(out2$WHEE.detected, out$detected)
    expect_identical(out2$WHEE.subset.sum.igg, out$subset.sum.igg)
    expect_identical(out2$WHEE.keep, out$keep)

    out3 <- quickAdtQc.se(se, subsets=list(igg=1:5, bar=6:10), flatten=FALSE)
    expect_identical(out$subset.sum.igg, out3$subset.sum$igg)
    expect_type(out3$subset.sum$bar, "double")

    out4 <- quickAdtQc.se(se, subsets=list(), meta.name=NULL)
    expect_null(metadata(out4)$qc)
})

test_that("quickAdtQc.se works with filtering", {
    set.seed(1399)
    mat <- matrix(rpois(10000, 1), ncol=100)
    se <- SummarizedExperiment(list(counts=mat))
    out <- quickAdtQc.se(se, subsets=list(), more.suggest.args=list(num.mads=0))
    expect_false(all(out$keep))

    filtered <- quickAdtQc.se(se, subsets=list(), filter.cells=TRUE, more.suggest.args=list(num.mads=0))
    expect_lt(ncol(filtered), ncol(se))
    expect_equal(ncol(filtered), sum(out$keep))
    expect_true(all(filtered$keep))
})

test_that("quickAdtQc.se works with character subsets", {
    ref <- quickAdtQc.se(se, subsets=list(Igg=1:10))

    # Check that we set the rownames when extracting the assay.
    copy <- se
    rownames(copy) <- sprintf("GENE_%s", seq_len(nrow(copy)))
    out <- quickAdtQc.se(copy, subsets=list(Igg=rownames(copy)[1:10]))

    expect_identical(colData(out), colData(ref))
})

test_that("quickAdtQc.se works with custom thresholds", {
    # Auto-fills the information that wasn't provided.
    fixed <- list(detected=pi, sum=1)
    out <- quickAdtQc.se(se, subsets=list(IgG=1:4), thresholds=fixed)

    expect_identical(metadata(out)$qc$thresholds$detected, pi)
    expect_null(metadata(out)$qc$thresholds$sum) # extra entry of thresholds is just ignored
    expect_identical(names(metadata(out)$qc$thresholds$subset.sum), "IgG")
    expect_identical(length(metadata(out)$qc$thresholds$subset.sum$IgG), 1L)

    # Used verbatim (some hacking required for renamed subsets).
    thresholds <- list(detected=100, subsets=list(IgG=5))
    out <- quickAdtQc.se(se, subsets=list(IgG=1:5), thresholds=I(thresholds))
    names(thresholds)[names(thresholds) == "subsets"] <- "subset.sum"
    expect_identical(metadata(out)$qc$thresholds, thresholds)
})

test_that("quickAdtQc.se overwrites existing entries", {
    copy <- se
    copy$sum <- "A"
    copy$detected <- "B"
    copy$subset.sum.Igg <- "C"

    copy <- quickAdtQc.se(copy, subsets=list(Igg=1:10))
    expect_type(copy$sum, "double")
    expect_type(copy$detected, "integer")
    expect_type(copy$subset.sum.Igg, "double")
})
