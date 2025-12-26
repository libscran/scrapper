# library(testthat); library(scrapper); source("test-se_quickAdtQc.R")

library(SummarizedExperiment)
mat <- matrix(rpois(1000, 1), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("quickAdtQc.se works as expected", {
    out <- quickAdtQc.se(se, subsets=list())
    expect_type(out$sum, "double")
    expect_type(out$detected, "integer")
    expect_type(out$keep, "logical")
    expect_type(metadata(out)$qc$thresholds, "list")

    out <- quickAdtQc.se(se, subsets=list(igg=1:5))
    expect_type(out$igg.sum, "double")

    out2 <- quickAdtQc.se(se, subsets=list(igg=1:5), output.prefix="WHEE.")
    expect_identical(out2$WHEE.sum, out$sum)
    expect_identical(out2$WHEE.detected, out$detected)
    expect_identical(out2$WHEE.igg.sum, out$igg.sum)
    expect_identical(out2$WHEE.keep, out$keep)

    out3 <- quickAdtQc.se(se, subsets=list(igg=1:5, bar=6:10), flatten=FALSE)
    expect_identical(out$igg.sum, out3$sum$igg)
    expect_type(out3$sum$bar, "double")

    out4 <- quickAdtQc.se(se, subsets=list(), meta.name=NULL)
    expect_null(metadata(out4)$qc)
})
