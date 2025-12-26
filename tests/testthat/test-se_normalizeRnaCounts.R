# library(testthat); library(scrapper); source("test-se_normalizeRnaCounts.R")

set.seed(12341234)
library(SummarizedExperiment)
mat <- matrix(rpois(1000, 5), ncol=10)
se <- SummarizedExperiment(list(counts=mat))

test_that("normalizeRnaCounts.se works as expected", {
    out <- normalizeRnaCounts.se(se)
    expect_type(out$sizeFactor, "double")
    expect_equal(mean(out$sizeFactor), 1)
    expect_s4_class(assay(out, "logcounts"), "DelayedArray")

    libs <- colSums(mat)
    out2 <- normalizeRnaCounts.se(se, size.factors=libs)
    expect_equal(out, out2)

    sf <- libs/mean(libs)
    out3 <- normalizeRnaCounts.se(se, size.factors=sf, center=FALSE)
    expect_equal(out, out3)

    out4 <- normalizeRnaCounts.se(se, factor.name=NULL)
    expect_null(out4$sizeFactor)
})
