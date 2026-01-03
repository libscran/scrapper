# library(testthat); library(scrapper); source("test-se_analyze.R")

set.seed(12345)
library(SingleCellExperiment)
mean <- runif(100, 0, 5)
mat <- matrix(rpois(10000, mean * 10), ncol=100)
se <- SummarizedExperiment(list(counts=mat))
rownames(se) <- paste0("GENE_", seq_len(nrow(se)))

# Just fiddling with some of the parameters to cut down the runtime.
default <- analyze.se(se, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)

test_that("analyze.se works correctly with a default run", {
    res <- default

    expect_identical(rownames(res$x), rownames(se))
    expect_true(all(res$x$keep))
    expect_type(res$x$keep, "logical")
    expect_null(res$x$combined.keep)

    expect_s4_class(counts(res$x), "DelayedMatrix")
    expect_s4_class(logcounts(res$x), "DelayedMatrix")

    expect_type(rowData(res$x)$residuals, "double")

    expect_true("PCA" %in% reducedDimNames(res$x))
    expect_true("TSNE" %in% reducedDimNames(res$x))
    expect_true("UMAP" %in% reducedDimNames(res$x))
    expect_false("combined" %in% reducedDimNames(res$x))
    expect_false("MNN" %in% reducedDimNames(res$x))

    expect_true(is.factor(res$x$graph.cluster))
    expect_identical(names(res$markers$rna), levels(res$x$graph.cluster))
})

test_that("analyze.se works correctly with no filtering", {
    res <- analyze.se(se, more.tsne.args=list(max.iterations=10), filter.cells=FALSE, more.umap.args=list(num.epochs=5), num.threads=1)
    expect_identical(dim(res$x), dim(se))
    expect_true(is.matrix(counts(res$x)))
})

test_that("analyze.se works correctly with k-means clustering", {
    res <- analyze.se(se, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), cluster.graph.output.name=NULL, kmeans.clusters=10, num.threads=1)
    expect_null(res$x$graph.cluster)
    expect_true(is.factor(res$x$kmeans.cluster))
    expect_identical(names(res$markers$rna), levels(res$x$kmeans.cluster))
})

test_that("analyze.se works correctly with no clustering", {
    res <- analyze.se(se, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), cluster.graph.output.name=NULL, num.threads=1)
    expect_null(res$x$graph.cluster)
    expect_null(res$x$kmeans.cluster)
    expect_null(res$markers)
})

test_that("analyze.se works correctly with blocking", {
    block <- rep(LETTERS[1:4], length.out=ncol(mat))
    res <- analyze.se(se, block=block, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)
    expect_true("MNN" %in% reducedDimNames(res$x))
    expect_type(res$x$block, "character")
    expect_identical(names(metadata(res$x)$qc$thresholds$sum), LETTERS[1:4])

    # Confirm that blocked coordinates are actually used downstream.
    downstream <- runAllNeighborSteps(as.matrix(t(reducedDim(res$x, "MNN"))), runTsne.args=list(max.iterations=10), runUmap.args=list(num.epochs=5), num.threads=1)
    expect_identical(downstream$clusterGraph$membership, res$x$graph.cluster)
    expect_identical(downstream$runTsne, reducedDim(res$x, "TSNE"))
})

test_that("analyze.se works correctly with combined RNA+ADT", {
    sce <- as(se, "SingleCellExperiment")
    amat <- matrix(rpois(1000, 10), ncol=100)
    rownames(amat) <- paste0("TAG_", seq_len(nrow(amat)))
    altExp(sce, "ADT") <- SummarizedExperiment(list(counts=amat))

    # Works with RNA + ADT. We use a different number of PCs for the ADTs to get some variety.
    res <- analyze.se(sce, adt.altexp="ADT", more.adt.pca.args=list(number=10), more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)

    expect_identical(rownames(altExp(res$x)), rownames(altExp(sce)))
    expect_type(altExp(res$x, "ADT")$keep, "logical")
    expect_type(res$x$combined.keep, "logical")
    expect_true(all(res$x$combined.keep))

    expect_s4_class(counts(altExp(res$x)), "DelayedMatrix")
    expect_s4_class(logcounts(altExp(res$x)), "DelayedMatrix")
    expect_null(rowData(altExp(res$x))$hvg)
    expect_identical(ncol(reducedDim(altExp(res$x), "PCA")), 10L)

    expect_type(rowData(res$x)$hvg, "logical")
    expect_identical(ncol(reducedDim(res$x, "PCA")), 25L)

    expect_identical(ncol(reducedDim(res$x, "combined")), 35L)
    expect_identical(names(metadata(res$x)$combined$main.scaling), "PCA")
    expect_identical(names(metadata(res$x)$combined$altexp.scaling$ADT), "PCA")

    expect_identical(names(res$markers$rna), levels(res$x$graph.cluster))
    expect_identical(names(res$markers$adt), levels(res$x$graph.cluster))

    # Also works with ADT + RNA where ADT is the main experiment.
    sce2 <- as(altExp(sce), "SingleCellExperiment")
    altExp(sce2, "RNA") <- se

    res2 <- analyze.se(sce2, rna.altexp="RNA", adt.altexp=NA, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)
    expect_null(rowData(res2$x)$hvg)
    expect_type(rowData(altExp(res2$x))$hvg, "logical")

    expect_identical(reducedDim(res2$x, "PCA"), reducedDim(altExp(res$x, "ADT"), "PCA"))
    expect_identical(reducedDim(altExp(res2$x, "RNA"), "PCA"), reducedDim(res$x, "PCA"))
    expect_identical(dim(reducedDim(res2$x, "combined")), dim(reducedDim(res$x, "combined")))

    expect_identical(res$markers, res2$markers)

    # Confirm that scaled coordinates are actually used downstream.
    ref <- runAllNeighborSteps(as.matrix(t(reducedDim(res$x, "combined"))), runTsne.args=list(max.iterations=10), runUmap.args=list(num.epochs=5), num.threads=1)
    expect_identical(ref$clusterGraph$membership, res$x$graph.cluster)
    expect_identical(ref$runTsne, reducedDim(res$x, "TSNE"))

    # Fails if neither ADT or RNA is supplied.
    expect_error(analyze.se(sce, use.rna.pcs=FALSE, use.adt.pcs=FALSE), "at least one")
})

test_that("analyze.se works correctly with ADT only", {
    amat <- matrix(rpois(1000, 10), ncol=100)
    rownames(amat) <- paste0("TAG_", seq_len(nrow(amat)))
    ase <- SummarizedExperiment(list(counts=amat))

    # Works with ADT in the main experiment.
    res <- analyze.se(ase, rna.altexp=NULL, adt.altexp=NA, more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)
    expect_type(res$x$keep, "logical")
    expect_true("PCA" %in% reducedDimNames(res$x))
    expect_false("combined" %in% reducedDimNames(res$x))
    expect_null(res$markers$rna)
    expect_identical(names(res$markers), "adt")

    # Works after ignoring RNA in the main experiment.
    # This checks that the use of a named character vector for 'target.embedding' is respected within analyze.se().
    sce <- as(se, "SingleCellExperiment")
    altExp(sce, "ADT") <- ase

    res2 <- analyze.se(sce, rna.altexp=NULL, adt.altexp="ADT", more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)
    expect_identical(reducedDim(altExp(res2$x, "ADT"), "PCA"), reducedDim(res$x, "PCA"))
    expect_identical(altExp(res2$x)$sizeFactor, res$x$sizeFactor)
    expect_identical(reducedDim(res2$x, "TSNE"), reducedDim(res$x, "TSNE")) # UMAP, TSNE, clusters still stored in the main experiment, though.
    expect_identical(reducedDim(res2$x, "UMAP"), reducedDim(res$x, "UMAP"))
    expect_identical(res2$x$graph.cluster, res$x$graph.cluster)
    expect_identical(res$markers, res2$markers)
})

test_that("analyze.se works correctly with the CRISPR modality", {
    sce <- as(se, "SingleCellExperiment")
    mat <- matrix(0, 20, 100)
    mat[cbind(sample(nrow(mat), ncol(mat), replace=TRUE), seq_len(ncol(mat)))] <- rpois(ncol(mat), 10)
    rownames(mat) <- paste0("GUIDE_", seq_len(nrow(mat)))
    altExp(sce, "CRISPR") <- SummarizedExperiment(list(counts=mat))

    res <- analyze.se(sce, crispr.altexp="CRISPR", more.tsne.args=list(max.iterations=10), more.umap.args=list(num.epochs=5), num.threads=1)
    expect_identical(rownames(altExp(res$x)), rownames(altExp(sce)))
    expect_type(altExp(res$x, "CRISPR")$keep, "logical")
    expect_type(res$x$combined.keep, "logical")
    expect_true(all(res$x$combined.keep))

    expect_s4_class(assay(altExp(res$x), "logcounts"), "DelayedMatrix")
    expect_identical(names(res$markers$crispr), levels(res$x$graph.cluster))
})

test_that(".extract_or_error works correctly", {
    se$foo <- runif(ncol(se))
    expect_identical(scrapper:::.extractOrError(colData(se), "foo"), se$foo)
    expect_error(scrapper:::.extractOrError(colData(se), "bar"), "no DataFrame")
})

test_that(".delayify_assays works correctly", {
    assay(se, "logcounts", withDimnames=FALSE) <- matrix(rpois(10000, mean * 10), ncol=100)
    sce <- as(se, "SingleCellExperiment")
    altExp(sce, "FOOBAR") <- se
    del <- scrapper:::.delayifyAssays(sce)

    expect_s4_class(assay(del, "counts"), "DelayedMatrix")
    expect_identical(as.matrix(assay(del, "counts")), assay(se, "counts"))

    expect_s4_class(assay(del, "logcounts"), "DelayedMatrix")
    expect_identical(as.matrix(assay(del, "logcounts")), assay(se, "logcounts"))

    expect_s4_class(assay(altExp(del), "counts"), "DelayedMatrix")
    expect_identical(as.matrix(assay(altExp(del), "counts")), assay(se, "counts"))
})

test_that(".define_single_target_embedding works as expected", {
    sce <- as(se, "SingleCellExperiment")
    altExp(sce, "FOOBAR") <- se
    expect_identical(scrapper:::.defineSingleTargetEmbedding(sce, NA, "PCA"), "PCA")
    expect_identical(scrapper:::.defineSingleTargetEmbedding(sce, 1, "PCA"), c(FOOBAR="PCA"))
    expect_identical(scrapper:::.defineSingleTargetEmbedding(sce, "STUFF", "PCA"), c(STUFF="PCA"))
})

test_that(".add_source_embedding_to_scale works as expected", {
    sce <- as(se, "SingleCellExperiment")
    altExp(sce, "FOOBAR") <- se

    everything <- list(main="A", altexp=list(YAY=2))
    expect_identical(scrapper:::.addSourceEmbeddingToScale(sce, NA, "PCA", everything), list(main=c("A", "PCA"), altexp=list(YAY=2)))
    expect_identical(scrapper:::.addSourceEmbeddingToScale(sce, 1, "PCA", everything), list(main="A", altexp=list(YAY=2, FOOBAR="PCA")))
    expect_identical(scrapper:::.addSourceEmbeddingToScale(sce, "STUFF", "PCA", everything), list(main="A", altexp=list(YAY=2, STUFF="PCA")))
})
