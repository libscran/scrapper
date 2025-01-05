# library(testthat); library(scrapper); source("test-analyze.R")

library(scRNAseq)
zeisel_sce <- fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays=TRUE)
zeisel_sce <- zeisel_sce[,1:500] # only using the first 500 cells for speed.

test_that("analyze works in the simple case", {
    # Lowering the MADs to check that filtering has some effect.
    res <- analyze(zeisel_sce, suggestRnaQcThresholds.args=list(num.mads=1), num.threads=2)

    expect_identical(res$combined.qc.filter, res$rna.qc.filter)
    expect_lt(sum(res$combined.qc.filter), ncol(zeisel_sce))
    expect_identical(nrow(res$rna.filtered), nrow(zeisel_sce))
    expect_lt(ncol(res$rna.filtered), ncol(zeisel_sce))

    expect_equal(mean(res$rna.size.factors), 1)
    expect_identical(dim(res$rna.filtered), dim(res$rna.normalized))

    expect_identical(ncol(res$rna.pca$components), ncol(res$rna.filtered))
    expect_identical(res$combined.pca, "rna.pca")
    expect_null(res$mnn.corrected)

    expect_identical(length(res$graph.clusters$membership), ncol(res$rna.normalized))
    expect_null(res$kmeans.clusters)
    nclusters <- length(unique(res$graph.clusters$membership))
    expect_gt(nclusters, 1)
    expect_identical(length(res$rna.markers$cohens.d), nclusters)

    expect_identical(ncol(res$tsne), 2L)
    expect_identical(ncol(res$umap), 2L)

    expect_null(res$adt.filtered)
    expect_null(res$adt.normalized)
    expect_null(res$adt.size.factors)
    expect_null(res$adt.pca)

    expect_null(res$crispr.filtered)
    expect_null(res$crispr.normalized)
    expect_null(res$crispr.size.factors)
    expect_null(res$crispr.pca)
})

test_that("analyze works with ADT data", {
    # We're going to pretend the spike-in data are ADTs, for testing purposes.
    # Again using a lower threshold to check that some filtering gets performed.
    res <- analyze(
        zeisel_sce,
        adt.x=altExp(zeisel_sce, "ERCC"),
        suggestRnaQcThresholds.args=list(num.mads=1),
        suggestAdtQcThresholds.args=list(num.mads=1),
        num.threads=2
    )

    expect_false(identical(res$combined.qc.filter, res$rna.qc.filter))
    expect_false(identical(res$combined.qc.filter, res$adt.qc.filter))
    expect_identical(res$combined.qc.filter, res$adt.qc.filter & res$rna.qc.filter)

    expect_identical(ncol(res$adt.filtered), ncol(res$rna.filtered))
    expect_identical(dim(res$adt.filtered), dim(res$adt.normalized))
    expect_equal(mean(res$adt.size.factors), 1)

    expect_identical(ncol(res$adt.pca$components), ncol(res$rna.pca$components))
    expect_identical(length(res$combined.pca$scaling), 2L)
    expect_identical(ncol(res$combined.pca$combined), ncol(res$adt.pca$components))
    expect_identical(nrow(res$combined.pca$combined), nrow(res$rna.pca$components) + nrow(res$adt.pca$components))

    expect_null(res$crispr.filtered)
    expect_null(res$crispr.normalized)
    expect_null(res$crispr.size.factors)
    expect_null(res$crispr.pca)
})

test_that("analyze works with CRISPR data", {
    # We're going to pretend the spike-in data are ADTs and CRISPR, for testing purposes.
    mock <- altExp(zeisel_sce, "ERCC")
    res <- analyze(rna.x=NULL, crispr.x=mock, adt.x=mock, num.threads=2)
    expect_identical(ncol(res$adt.filtered), ncol(res$crispr.filtered))
    expect_identical(dim(res$crispr.filtered), dim(res$crispr.normalized))
    expect_equal(mean(res$crispr.size.factors), 1)

    expect_identical(ncol(res$adt.pca$components), ncol(res$crispr.pca$components))
    expect_identical(length(res$combined.pca$scaling), 2L)
    expect_identical(ncol(res$combined.pca$combined), ncol(res$adt.pca$components))
    expect_identical(nrow(res$combined.pca$combined), nrow(res$crispr.pca$components) + nrow(res$adt.pca$components))

    # Checking that this all works even if no RNA data is supplied$
    expect_null(res$rna.filtered)
    expect_null(res$rna.normalized)
    expect_null(res$rna.size.factors)
    expect_null(res$rna.pca)
})

test_that("analyze works with blocking", {
    # Using the tissue as the blocking factor for testing purposes.
    block <- zeisel_sce$tissue
    levels <- sort(unique(block))
    res <- analyze(rna.x=assay(zeisel_sce), block=block, num.threads=2)

    expect_true(!is.null(res$mnn.corrected))
    expect_identical(sort(res$mnn.corrected$merge.order), levels)
    expect_identical(length(res$rna.qc.thresholds$sum), length(levels))
    expect_identical(nrow(res$rna.pca$center), length(levels))
    expect_identical(length(res$rna.gene.variances$per.block), length(levels))
    expect_false(isTRUE(all.equal(mean(res$rna.size.factors), 1))) # as the size factors are only scaled to a mean of 1 in the lowest-coverage block.
})

test_that("analyze works without filtering", {
    # Lowering the MADs to check that filtering would some effect, if not for filter.cells=FALSE.
    res <- analyze(
        assay(zeisel_sce),
        filter.cells=FALSE,
        suggestRnaQcThresholds.args=list(num.mads=1),
        num.threads=2
    )
    expect_lt(sum(res$combined.qc.filter), ncol(zeisel_sce))
    expect_identical(dim(res$rna.filtered), dim(zeisel_sce))
})

test_that("analyze works with k-means", {
    res <- analyze(
        assay(zeisel_sce),
        kmeans.clusters=13,
        num.threads=2
    )
    expect_identical(length(unique(res$kmeans.clusters$clusters)), 13L)
    expect_identical(length(res$rna.markers$auc), length(unique(res$graph.clusters$membership)))

    res <- analyze(
        assay(zeisel_sce), 
        kmeans.clusters=13,
        clusters.for.markers="kmeans",
        num.threads=2
    )
    expect_identical(length(res$rna.markers$auc), 13L)
})

test_that("analyze works without any clustering", {
    res <- analyze(
        assay(zeisel_sce),
        clusterGraph.args=NULL,
        num.threads=2
    )

    # No clusters, no markers.
    expect_null(res$kmeans.clusters)
    expect_null(res$graph.clusters)
    expect_null(res$rna.markers)
})
