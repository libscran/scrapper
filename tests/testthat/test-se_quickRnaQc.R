# library(testthat); library(scrapper); source("test-se_quickRnaQc.R")

library(SingleCellExperiment)
mat <- matrix(rpois(1000, 1), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("quickRnaQc.se works as expected", {
    out <- quickRnaQc.se(se, subsets=list())
    expect_identical(ncol(out), ncol(se))
    expect_type(out$sum, "double")
    expect_type(out$detected, "integer")
    expect_type(out$keep, "logical")
    expect_type(metadata(out)$qc$thresholds, "list")

    out <- quickRnaQc.se(se, subsets=list(mt=1:5))
    expect_type(out$subset.proportion.mt, "double")

    out2 <- quickRnaQc.se(se, subsets=list(mt=1:5), output.prefix="WHEE.")
    expect_identical(out2$WHEE.sum, out$sum)
    expect_identical(out2$WHEE.detected, out$detected)
    expect_identical(out2$WHEE.subset.proportion.mt, out$subset.proportion.mt)
    expect_identical(out2$WHEE.keep, out$keep)

    out3 <- quickRnaQc.se(se, subsets=list(mt=1:5, bar=6:10), flatten=FALSE)
    expect_identical(out$subset.proportion.mt, out3$subset.proportion$mt)
    expect_type(out3$subset.proportion$bar, "double")

    out4 <- quickRnaQc.se(se, subsets=list(), meta.name=NULL)
    expect_null(metadata(out4)$qc)
})

test_that("quickRnaQc.se works with filtering", {
    set.seed(102938)
    mat <- matrix(rpois(10000, 1), ncol=100)
    se <- SummarizedExperiment(list(counts=mat))
    out <- quickRnaQc.se(se, subsets=list(), more.suggest.args=list(num.mads=0))
    expect_false(all(out$keep))

    filtered <- quickRnaQc.se(se, subsets=list(), filter.cells=TRUE, more.suggest.args=list(num.mads=0))
    expect_lt(ncol(filtered), ncol(se))
    expect_equal(ncol(filtered), sum(out$keep))
    expect_true(all(filtered$keep))
})

test_that("quickRnaQc.se works with character subsets", {
    ref <- quickRnaQc.se(se, subsets=list(Mt=1:10))

    # Check that we set the rownames when extracting the assay.
    copy <- se
    rownames(copy) <- sprintf("GENE_%s", seq_len(nrow(copy)))
    out <- quickRnaQc.se(copy, subsets=list(Mt=rownames(copy)[1:10]))

    expect_identical(colData(out), colData(ref))
})

test_that("quickRnaQc.se works with alternative experiments", {
    sce <- as(se, "SingleCellExperiment")
    altExp(sce, "ERCC") <- se[1:20,]

    out <- quickRnaQc.se(sce, subsets=list(), altexp.proportions="ERCC")
    expect_type(out$subset.proportion.ERCC, "double")
    expect_equal(1 / (1 + out$sum / altExp(out)$sum), out$subset.proportion.ERCC)

    out <- quickRnaQc.se(sce, subsets=list(), altexp.proportions="ERCC", output.prefix="WHEE.")
    expect_type(out$WHEE.subset.proportion.ERCC, "double")
    expect_equal(1 / (1 + out$WHEE.sum / altExp(out, "ERCC")$WHEE.sum), out$WHEE.subset.proportion.ERCC)
})

test_that("quickRnaQc.se works with custom thresholds", {
    # Auto-fills the information that wasn't provided.
    fixed <- list(sum=pi)
    out <- quickRnaQc.se(se, subsets=list(Mito=1:5), thresholds=fixed)

    expect_identical(metadata(out)$qc$thresholds$sum, pi)
    expect_identical(length(metadata(out)$qc$thresholds$detected), 1L)
    expect_identical(names(metadata(out)$qc$thresholds$subset.proportion), "Mito")
    expect_identical(length(metadata(out)$qc$thresholds$subset.proportion$Mito), 1L)

    # Used verbatim (some hacking required for renamed subsets).
    thresholds <- list(sum=pi, detected=100, subsets=list(Mito=5))
    out <- quickRnaQc.se(se, subsets=list(Mito=1:5), thresholds=I(thresholds))
    names(thresholds)[names(thresholds) == "subsets"] <- "subset.proportion"
    expect_identical(metadata(out)$qc$thresholds, thresholds)
})

test_that(".substituteThresholds works correctly", {
    expect_identical(scrapper:::.substituteThresholds(list(a=1, b=2, d=5), list(a=2, b=3, c=4)), list(a=2, b=3, d=5))
    expect_identical(scrapper:::.substituteThresholds(list(a=1, b=c(x=2, y=4), d=5), list(a=2, b=3, c=4)), list(a=2, b=c(x=3, y=3), d=5))
    expect_identical(scrapper:::.substituteThresholds(list(a=1, b=c(x=2, y=4), d=5), list(a=2, b=c(x=3, y=6), c=4)), list(a=2, b=c(x=3, y=6), d=5))
    expect_error(scrapper:::.substituteThresholds(list(a=1, b=2, d=5), list(a=2, b=c(x=1), c=4)), "unexpected names")
    expect_error(scrapper:::.substituteThresholds(list(a=1, b=2, d=5), list(a=2, b=1:2, c=4)), "expected a scalar")
})

test_that("quickRnaQc.se overwrites existing entries", {
    copy <- se
    copy$sum <- "A"
    copy$detected <- "B"
    copy$subset.proportion.Mt <- "C"

    copy <- quickRnaQc.se(copy, subsets=list(Mt=1:10))
    expect_type(copy$sum, "double")
    expect_type(copy$detected, "integer")
    expect_type(copy$subset.proportion.Mt, "double")
})
