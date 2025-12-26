# library(testthat); library(scrapper); source("test-se_normalizeAdtCounts.R")

set.seed(12341234)
library(SummarizedExperiment)
mat <- matrix(rpois(500, 50), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("normalizeAdtCounts.se works as expected", {
    out <- normalizeAdtCounts.se(se)
    expect_type(out$sizeFactor, "double")
    expect_equal(mean(out$sizeFactor), 1)
    expect_s4_class(assay(out, "logcounts"), "DelayedArray")

    raw.sf <- computeClrm1Factors(mat)
    out2 <- normalizeAdtCounts.se(se, size.factors=raw.sf)
    expect_equal(out, out2)

    sf <- raw.sf/mean(raw.sf)
    out3 <- normalizeAdtCounts.se(se, size.factors=sf, center=FALSE)
    expect_equal(out, out3)

    out4 <- normalizeAdtCounts.se(se, factor.name=NULL)
    expect_null(out4$sizeFactor)
})
