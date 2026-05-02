# library(testthat); library(scrapper); source("test-se_normalizeRnaCountsWithSpikeIns.R")

set.seed(96969696)
library(SingleCellExperiment)
mat <- matrix(rpois(1000, 5), ncol=10)
sce <- SingleCellExperiment(list(counts=mat))
emat <- matrix(rpois(200, 5), ncol=10)
altExp(sce, "ERCC") <- SingleCellExperiment(list(counts=emat))
smat <- matrix(rpois(500, 5), ncol=10)
altExp(sce, "SIRV") <- SingleCellExperiment(list(counts=smat))

test_that("normalizeRnaCountsWithSpikeIns works as expected", {
    e.centered <- centerSizeFactors(colSums(mat)) 
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"))

    expect_equal(out$sizeFactor, e.centered)
    expect_equal(altExp(out, "ERCC")$sizeFactor, centerSizeFactors(colSums(emat)))
    expect_equal(altExp(out, "SIRV")$sizeFactor, centerSizeFactors(colSums(smat)))

    expect_s4_class(assay(out, "logcounts"), "DelayedArray")
    expect_s4_class(assay(altExp(out, "ERCC"), "logcounts"), "DelayedArray")
    expect_s4_class(assay(altExp(out, "SIRV"), "logcounts"), "DelayedArray")

    # Works with nothing.
    none <- normalizeRnaCountsWithSpikeIns.se(sce, NULL)
    expect_equal(none$sizeFactor, e.centered)
    expect_null(altExp(none, "ERCC")$sizeFactor)
    expect_null(altExp(none, "SIRV")$sizeFactor)

    # Works with other ways of specifying the alternative experiments.
    alt <- normalizeRnaCountsWithSpikeIns.se(sce, c(SIRV = 1, ERCC = 1))
    expect_identical(out, alt)
})

test_that("normalizeRnaCountsWithSpikeIns works with manually supplied factors", {
    sf <- runif(10)
    esf <- runif(10)
    ssf <- runif(10)
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"), endogenous.factors=sf, spike.factors=list(ERCC=esf, SIRV=ssf))

    expect_equal(out$sizeFactor, centerSizeFactors(sf))
    expect_equal(altExp(out, "ERCC")$sizeFactor, centerSizeFactors(esf))
    expect_equal(altExp(out, "SIRV")$sizeFactor, centerSizeFactors(ssf))
})

test_that("normalizeRnaCountsWithSpikeIns can disable storage of the size factors", {
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"), factor.name=NULL)
    expect_null(out$sizeFactor)
    expect_null(altExp(out, "ERCC")$sizeFactor)
    expect_null(altExp(out, "SIRV")$sizeFactor)
})

test_that("normalizeRnaCountsWithSpikeIns supports normalization on spike-ins", {
    e.centered <- centerSizeFactors(colSums(emat))
    s.centered <- centerSizeFactors(colSums(smat))

    out <- normalizeRnaCountsWithSpikeIns.se(sce, "ERCC", use.spike.ins.for.endogenous=TRUE)
    expect_equal(out$sizeFactor, e.centered)
    expect_equal(altExp(out, "ERCC")$sizeFactor, e.centered)

    # Works with multiple entries.
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"), use.spike.ins.for.endogenous=TRUE)
    expect_equal(out$sizeFactor, e.centered)
    expect_equal(altExp(out, "ERCC")$sizeFactor, e.centered)
    expect_equal(altExp(out, "SIRV")$sizeFactor, s.centered)

    # Works with a named altexp.
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"), use.spike.ins.for.endogenous="SIRV")
    expect_equal(out$sizeFactor, s.centered)
    expect_equal(altExp(out, "ERCC")$sizeFactor, e.centered)
    expect_equal(altExp(out, "SIRV")$sizeFactor, s.centered)
})

test_that("normalizeRnaCountsWithSpikeIns works with blocking", {
    b <- rep(1:3, length.out=ncol(sce))
    out <- normalizeRnaCountsWithSpikeIns.se(sce, c("ERCC", "SIRV"), block=b)

    sf.means <- vapply(split(out$sizeFactor, b), mean, 0)
    expect_equal(min(sf.means), 1)
    esf.means <- vapply(split(altExp(out, "ERCC")$sizeFactor, b), mean, 0)
    expect_equal(sf.means, esf.means)
    ssf.means <- vapply(split(altExp(out, "SIRV")$sizeFactor, b), mean, 0)
    expect_equal(sf.means, ssf.means)
})
