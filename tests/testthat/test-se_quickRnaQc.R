# library(testthat); library(scrapper); source("test-se_quickRnaQc.R")

library(SingleCellExperiment)
mat <- matrix(rpois(1000, 1), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("quickRnaQc.se works as expected", {
    out <- quickRnaQc.se(se, subsets=list())
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
