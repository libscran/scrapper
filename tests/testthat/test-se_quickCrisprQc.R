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
