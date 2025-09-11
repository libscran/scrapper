#' Analyze single-cell data
#'
#' Execute a simple single-cell analysis pipeline, starting from a count matrix and ending with clusters, visualizations and markers.
#' This also supports integration of multiple modalities and correction of batch effects.
#'
#' @param rna.x Matrix-like object containing RNA counts.
#' This should have the same number of columns as the other \code{*.x} arguments.
#'
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} instance containing such a matrix in its \code{rna.assay}.
#'
#' Alternatively \code{NULL}, if no RNA counts are available.
#' @param adt.x Matrix-like object containing ADT counts.
#' This should have the same number of columns as the other \code{*.x} arguments.
#'
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} instance containing such a matrix in its \code{adt.assay}.
#'
#' Alternatively \code{NULL}, if no ADT counts are available.
#' @param crispr.x Matrix-like object containing ADT counts.
#' This should have the same number of columns as the other \code{*.x} arguments.
#'
#' Alternatively, a \link[SummarizedExperiment]{SummarizedExperiment} instance containing such a matrix in its \code{crispr.assay}.
#'
#' Alternatively \code{NULL}, if no ADT counts are available.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in the \code{*_x} matrices.
#' Alternatively \code{NULL}, if all cells are from the same block.
#' @param rna.subsets Gene subsets for quality control, typically used for mitochondrial genes.
#' See the \code{subsets} arguments in \code{\link{computeRnaQcMetrics}} for details.
#' @param adt.subsets ADT subsets for quality control, typically used for IgG controls.
#' See the \code{subsets} arguments in \code{\link{computeAdtQcMetrics}} for details.
#' @param suggestRnaQcThresholds.args Named list of arguments to pass to \code{\link{suggestRnaQcThresholds}}.
#' @param suggestAdtQcThresholds.args Named list of arguments to pass to \code{\link{suggestAdtQcThresholds}}.
#' @param suggestCrisprQcThresholds.args Named list of arguments to pass to \code{\link{suggestCrisprQcThresholds}}.
#' @param filter.cells Logical scalar indicating whether to filter the count matrices to only retain high-quality cells in all modalities.
#' If \code{FALSE}, QC metrics and thresholds are still computed but are not used to filter the count matrices.
#' @param centerSizeFactors.args Named list of arguments to pass to \code{\link{centerSizeFactors}}.
#' @param computeClrm1Factors.args Named list of arguments to pass to \code{\link{computeClrm1Factors}}.
#' Only used if \code{adt.x} is provided.
#' @param normalizeCounts.args Named list of arguments to pass to \code{\link{normalizeCounts}}.
#' @param modelGeneVariances.args Named list of arguments to pass to \code{\link{modelGeneVariances}}.
#' Only used if \code{rna.x} is provided.
#' @param chooseHighlyVariableGenes.args Named list of arguments to pass to \code{\link{chooseHighlyVariableGenes}}.
#' Only used if \code{rna.x} is provided.
#' @param runPca.args Named list of arguments to pass to \code{\link{runPca}}.
#' @param use.rna.pcs Logical scalar indicating whether to use the RNA-derived PCs for downstream steps (i.e., clustering, visualization).
#' Only used if \code{rna.x} is provided.
#' @param use.adt.pcs Logical scalar indicating whether to use the ADT-derived PCs for downstream steps (i.e., clustering, visualization).
#' Only used if \code{adt.x} is provided.
#' @param use.crispr.pcs Logical scalar indicating whether to use the CRISPR-derived PCs for downstream steps (i.e., clustering, visualization).
#' Only used if \code{crispr.x} is provided.
#' @param scaleByNeighbors.args Named list of arguments to pass to \code{\link{scaleByNeighbors}}.
#' Only used if multiple modalities are available and their corresponding \code{use.*.pcs} arguments are \code{TRUE}.
#' @param correctMnn.args Named list of arguments to pass to \code{\link{correctMnn}}.
#' Only used if \code{block} is supplied.
#' @param runTsne.args Named list of arguments to pass to \code{\link{runTsne}}.
#' If \code{NULL}, t-SNE is not performed.
#' @param runUmap.args Named list of arguments to pass to \code{\link{runUmap}}.
#' If \code{NULL}, UMAP is not performed.
#' @param buildSnnGraph.args Named list of arguments to pass to \code{\link{buildSnnGraph}}.
#' Ignored if \code{clusterGraph.args = NULL}.
#' @param clusterGraph.args Named list of arguments to pass to \code{\link{clusterGraph}}.
#' If \code{NULL}, graph-based clustering is not performed.
#' @param runAllNeighborSteps.args Named list of arguments to pass to \code{\link{runAllNeighborSteps}}.
#' @param kmeans.clusters Integer scalar specifying the number of clusters to use in k-means clustering.
#' If \code{NULL}, k-means clustering is not performed.
#' @param clusterKmeans.args Named list of arguments to pass to \code{\link{clusterKmeans}}.
#' Ignored if \code{kmeans.clusters = NULL}.
#' @param clusters.for.markers Character vector of clustering algorithms (either \code{"graph"} or \code{"kmeans"}, specifying the clustering to be used for marker detection.
#' The first available clustering will be chosen.
#' @param scoreMarkers.args Named list of arguments to pass to \code{\link{scoreMarkers}}.
#' Ignored if no suitable clusterings are available.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} instance specifying the nearest-neighbor search algorithm to use.
#' @param rna.assay Integer scalar or string specifying the assay to use if \code{rna.x} is a \link[SummarizedExperiment]{SummarizedExperiment}.
#' @param adt.assay Integer scalar or string specifying the assay to use if \code{adt.x} is a \link[SummarizedExperiment]{SummarizedExperiment}.
#' @param crispr.assay Integer scalar or string specifying the assay to use if \code{crispr.x} is a \link[SummarizedExperiment]{SummarizedExperiment}.
#' @param num.threads Integer scalar specifying the number of threads to use in each step.
#'
#' @return List containing the results of the entire analysis:
#' \describe{
#' \item{\code{rna.qc.metrics}:}{Results of \code{\link{computeRnaQcMetrics}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.qc.thresholds}:}{Results of \code{\link{suggestRnaQcThresholds}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.qc.filter}:}{Results of \code{\link{filterRnaQcMetrics}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.qc.metrics}:}{Results of \code{\link{computeAdtQcMetrics}}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.qc.thresholds}:}{Results of \code{\link{suggestAdtQcThresholds}}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.qc.filter}:}{Results of \code{\link{filterAdtQcMetrics}}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.qc.metrics}:}{Results of \code{\link{computeCrisprQcMetrics}}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.qc.thresholds}:}{Results of \code{\link{suggestCrisprQcThresholds}}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.qc.filter}:}{Results of \code{\link{filterCrisprQcMetrics}}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{combined.qc.filter}:}{Logical vector indicating which cells are of high quality and should be retained for downstream analyses.}
#' \item{\code{rna.filtered}:}{Matrix of RNA counts that has been filtered to only contain the high-quality cells in \code{combined.qc.filter}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.filtered}:}{Matrix of ADT counts that has been filtered to only contain the high-quality cells in \code{combined.qc.filter}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.filtered}:}{Matrix of CRISPR counts that has been filtered to only contain the high-quality cells in \code{combined.qc.filter}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.size.factors}:}{Size factors for the RNA count matrix, derived from the sum of counts for each cell and centered with \code{\link{centerSizeFactors}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.normalized}:}{Matrix of (log-)normalized expression values derived from RNA counts, as computed by \code{\link{normalizeCounts}} using \code{rna.size.factors}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.size.factors}:}{Size factors for the ADT count matrix, computed by \code{\link{computeClrm1Factors}} and centered with \code{\link{centerSizeFactors}}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.normalized}:}{Matrix of (log-)normalized expression values derived from ADT counts, as computed by \code{\link{normalizeCounts}} using \code{adt.size.factors}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.size.factors}:}{Size factors for the CRISPR count matrix, derived from the sum of counts for each cell and centered with \code{\link{centerSizeFactors}}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.normalized}:}{Matrix of (log-)normalized expression values derived from CRISPR counts, as computed by \code{\link{normalizeCounts}} using \code{crispr.size.factors}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.gene.variances}:}{Results of \code{\link{modelGeneVariances}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.highly.variable.genes}:}{Results of \code{\link{chooseHighlyVariableGenes}}.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{rna.pca}:}{Results of calling \code{\link{runPca}} on \code{rna.normalized} with the \code{rna.highly.variable.genes} subset.
#' If RNA data is not available, this is set to \code{NULL} instead.}
#' \item{\code{adt.pca}:}{Results of calling \code{\link{runPca}} on \code{adt.normalized}.
#' If ADT data is not available, this is set to \code{NULL} instead.}
#' \item{\code{crispr.pca}:}{Results of calling \code{\link{runPca}} on \code{crispr.normalized}.
#' If CRISPR data is not available, this is set to \code{NULL} instead.}
#' \item{\code{combined.pca}:}{If only one modality is used for the downstream analysis, this is a string specifying the list element containing the components to be used, e.g., \code{"rna.pca"}.
#' If multiple modalities are to be combined for downstream analysis, this contains the results of \code{\link{scaleByNeighbors}} on the PCs of those modalities.}
#' \item{\code{block}:}{Vector or factor containing the blocking factor for all cells (after filtering, if \code{filter.cells = TRUE}).
#' This is set to \code{NULL} if no blocking factor was supplied.}
#' \item{\code{mnn.corrected}:}{Results of \code{\link{correctMnn}} on the PCs in or referenced by \code{combined.pca}.
#' If no blocking factor is supplied, this is set to \code{NULL} instead.}
#' \item{\code{tsne}:}{Results of \code{\link{runTsne}}.
#' This is \code{NULL} if t-SNE was not performed.}
#' \item{\code{umap}:}{Results of \code{\link{runUmap}}.
#' This is \code{NULL} if UMAP was not performed.}
#' \item{\code{snn.graph}:}{Results of \code{\link{buildSnnGraph}}.
#' This is \code{NULL} if graph-based clustering was not performed, or if \code{return.graph=FALSE} in \code{\link{runAllNeighborSteps}}.}
#' \item{\code{graph.clusters}:}{Results of \code{\link{clusterGraph}}.
#' This is \code{NULL} if graph-based clustering was not performed.}
#' \item{\code{kmeans.clusters}:}{Results of \code{\link{clusterKmeans}}.
#' This is \code{NULL} if k-means clustering was not performed.}
#' \item{\code{clusters}:}{Integer vector containing the cluster assignment for each cell (after filtering, if \code{filter.cells = TRUE}).
#' This may be derived from \code{graph.clusters} or \code{kmeans.clusters} depending on the choice of \code{clusters.for.markers}.
#' If no suitable clusterings are available, this is set to NULL.}
#' \item{\code{rna.markers}:}{Results of calling \code{\link{scoreMarkers}} on \code{rna.normalized}.
#' This is \code{NULL} if RNA data is not available or no suitable clusterings are available.}
#' \item{\code{adt.markers}:}{Results of calling \code{\link{scoreMarkers}} on \code{adt.normalized}.
#' This is \code{NULL} if ADT data is not available or no suitable clusterings are available.}
#' \item{\code{crispr.markers}:}{Results of calling \code{\link{scoreMarkers}} on \code{crispr.normalized}.
#' This is \code{NULL} if CRISPR data is not available or no suitable clusterings are available.}
#' }
#'
#' @seealso
#' \code{\link{convertAnalyzeResults}}, to convert the results into a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @author Aaron Lun
#' @examples
#' library(scRNAseq)
#' sce <- fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays=TRUE)
#' sce <- sce[,1:500] # smaller dataset for a faster runtime for R CMD check. 
#' res <- analyze(
#'     sce, 
#'     rna.subsets=list(mito=grep("^mt-", rownames(sce))),
#'     num.threads=2 # keep R CMD check happy
#' )
#' str(res)
#' convertAnalyzeResults(res)
#'
#' @export
#' @importFrom methods is
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom DelayedArray DelayedArray
analyze <- function(
    rna.x,
    adt.x = NULL,
    crispr.x = NULL,
    block = NULL,
    rna.subsets = list(),
    adt.subsets = list(),
    suggestRnaQcThresholds.args = list(),
    suggestAdtQcThresholds.args = list(),
    suggestCrisprQcThresholds.args = list(),
    filter.cells = TRUE,
    centerSizeFactors.args = list(),
    computeClrm1Factors.args = list(),
    normalizeCounts.args = list(),
    modelGeneVariances.args = list(),
    chooseHighlyVariableGenes.args = list(),
    runPca.args = list(),
    use.rna.pcs = TRUE,
    use.adt.pcs = TRUE,
    use.crispr.pcs = TRUE,
    scaleByNeighbors.args = list(),
    correctMnn.args = list(),
    runUmap.args = list(),
    runTsne.args = list(),
    buildSnnGraph.args = list(),
    clusterGraph.args = list(),
    runAllNeighborSteps.args = list(),
    kmeans.clusters = NULL,
    clusterKmeans.args = list(),
    clusters.for.markers = c("graph", "kmeans"),
    scoreMarkers.args = list(),
    BNPARAM = AnnoyParam(),
    rna.assay = 1L,
    adt.assay = 1L,
    crispr.assay = 1L,
    num.threads = 3L)
{
    store <- list()

    try.se.extract <- function(x, assay) {
        if (is(x, "SummarizedExperiment")) {
            SummarizedExperiment::assay(x, assay)
        } else {
            x
        }
    }

    all.ncols <- integer()
    if (!is.null(rna.x)) {
        all.ncols <- c(all.ncols, ncol(rna.x))
        rna.x <- try.se.extract(rna.x, rna.assay)
    }
    if (!is.null(adt.x)) {
        all.ncols <- c(all.ncols, ncol(adt.x))
        adt.x <- try.se.extract(adt.x, adt.assay)
    }
    if (!is.null(crispr.x)) {
        all.ncols <- c(all.ncols, ncol(crispr.x))
        crispr.x <- try.se.extract(crispr.x, crispr.assay)
    }

    ncells <- unique(all.ncols)
    if (length(ncells) > 1) {
        stop("'*.x' have differing numbers of columns")
    }
    if (length(ncells) == 0) {
        stop("at least one '*.x' must be non-NULL")
    }

    ############ Quality control o(*°▽°*)o #############

    if (!is.null(rna.x)) {
        store[["rna.qc.metrics"]] <- computeRnaQcMetrics(rna.x, subsets=rna.subsets, num.threads=num.threads)
        store[["rna.qc.thresholds"]] <- do.call(suggestRnaQcThresholds, c(list(store[["rna.qc.metrics"]], block=block), suggestRnaQcThresholds.args))
        store[["rna.qc.filter"]] <- filterRnaQcMetrics(store[["rna.qc.thresholds"]], store[["rna.qc.metrics"]], block=block)
    } else {
        store["rna.qc.metrics"] <- list(NULL)
        store["rna.qc.thresholds"] <- list(NULL)
        store["rna.qc.filter"] <- list(NULL)
    }

    if (!is.null(adt.x)) {
        store[["adt.qc.metrics"]] <- computeAdtQcMetrics(adt.x, subsets=adt.subsets, num.threads=num.threads)
        store[["adt.qc.thresholds"]] <- do.call(suggestAdtQcThresholds, c(list(store[["adt.qc.metrics"]], block=block), suggestAdtQcThresholds.args))
        store[["adt.qc.filter"]] <- filterAdtQcMetrics(store[["adt.qc.thresholds"]], store[["adt.qc.metrics"]], block=block)
    } else {
        store["adt.qc.metrics"] <- list(NULL)
        store["adt.qc.thresholds"] <- list(NULL)
        store["adt.qc.filter"] <- list(NULL)
    }

    if (!is.null(crispr.x)) {
        store[["crispr.qc.metrics"]] <- computeCrisprQcMetrics(crispr.x, num.threads=num.threads)
        store[["crispr.qc.thresholds"]] <- do.call(suggestCrisprQcThresholds, c(list(store[["crispr.qc.metrics"]], block=block), suggestCrisprQcThresholds.args))
        store[["crispr.qc.filter"]] <- filterCrisprQcMetrics(store[["crispr.qc.thresholds"]], store[["crispr.qc.metrics"]], block=block)
    } else {
        store["crispr.qc.metrics"] <- list(NULL)
        store["crispr.qc.thresholds"] <- list(NULL)
        store["crispr.qc.filter"] <- list(NULL)
    }

    # Combining all filters.
    combined.qc.filter <- !logical(ncells)
    for (mod in c("rna", "adt", "crispr")) {
        keep <- store[[paste0(mod, ".qc.filter")]]
        if (!is.null(keep)) {
            combined.qc.filter <- combined.qc.filter & keep
        }
    }
    store[["combined.qc.filter"]] <- combined.qc.filter

    if (!is.null(rna.x)) {
        store[["rna.filtered"]] <- DelayedArray(rna.x)
    } else {
        store["rna.filtered"] <- list(NULL)
    }
    if (!is.null(adt.x)) {
        store[["adt.filtered"]] <- DelayedArray(adt.x)
    } else {
        store["adt.filtered"] <- list(NULL)
    }
    if (!is.null(crispr.x)) {
        store[["crispr.filtered"]] <- DelayedArray(crispr.x)
    } else {
        store["crispr.filtered"] <- list(NULL)
    }

    if (filter.cells) {
        if (!is.null(rna.x)) {
            store[["rna.filtered"]] <- store[["rna.filtered"]][,combined.qc.filter,drop=FALSE]
        }
        if (!is.null(adt.x)) {
            store[["adt.filtered"]] <- store[["adt.filtered"]][,combined.qc.filter,drop=FALSE]
        }
        if (!is.null(crispr.x)) {
            store[["crispr.filtered"]] <- store[["crispr.filtered"]][,combined.qc.filter,drop=FALSE]
        }
        if (!is.null(block)) {
            block <- block[combined.qc.filter]
        }
    }

    ############ Normalization ( ꈍᴗꈍ) #############

    maybe.filter <- function(x) {
        if (filter.cells && !all(combined.qc.filter)) {
            x[combined.qc.filter]
        } else {
            x
        }
    }

    if (!is.null(rna.x)) {
        store[["rna.size.factors"]] <- do.call(centerSizeFactors, c(list(maybe.filter(store[["rna.qc.metrics"]]$sum), block=block), centerSizeFactors.args))
        store[["rna.normalized"]] <- do.call(normalizeCounts, c(list(store[["rna.filtered"]], store[["rna.size.factors"]]), normalizeCounts.args))
    } else {
        store["rna.size.factors"] <- list(NULL)
        store["rna.normalized"] <- list(NULL)
    }

    if (!is.null(adt.x)) {
        raw.adt.sf <- do.call(computeClrm1Factors, c(list(store[["adt.filtered"]], num.threads=num.threads), computeClrm1Factors.args))
        store[["adt.size.factors"]] <- do.call(centerSizeFactors, c(list(raw.adt.sf, block=block), centerSizeFactors.args))
        store[["adt.normalized"]] <- do.call(normalizeCounts, c(list(store[["adt.filtered"]], store[["adt.size.factors"]]), normalizeCounts.args))
    } else {
        store["adt.size.factors"] <- list(NULL)
        store["adt.normalized"] <- list(NULL)
    }

    if (!is.null(crispr.x)) {
        store[["crispr.size.factors"]] <- do.call(centerSizeFactors, c(list(maybe.filter(store[["crispr.qc.metrics"]]$sum), block=block), centerSizeFactors.args))
        store[["crispr.normalized"]] <- do.call(normalizeCounts, c(list(store[["crispr.filtered"]], store[["crispr.size.factors"]]), normalizeCounts.args))
    } else {
        store["crispr.size.factors"] <- list(NULL)
        store["crispr.normalized"] <- list(NULL)
    }

    ############ Variance modelling (～￣▽￣)～ #############

    if (!is.null(rna.x)) {
        store[["rna.gene.variances"]] <- do.call(modelGeneVariances, c(list(store[["rna.normalized"]], block=block, num.threads=num.threads), modelGeneVariances.args))
        store[["rna.highly.variable.genes"]] <- do.call(chooseHighlyVariableGenes, c(list(store[["rna.gene.variances"]]$statistics$residuals), chooseHighlyVariableGenes.args))
    } else {
        store["rna.gene.variances"] <- list(NULL)
        store["rna.highly.variable.genes"] <- list(NULL)
    }

    ############ Principal components analysis \(>⩊<)/ #############

    if (!is.null(rna.x)) {
        hvg.mat <- store[["rna.normalized"]][store[["rna.highly.variable.genes"]],,drop=FALSE]
        store[["rna.pca"]] <- do.call(runPca, c(list(hvg.mat, block=block, num.threads=num.threads), runPca.args))
    } else {
        store["rna.pca"] <- list(NULL)
    }

    if (!is.null(adt.x)) {
        store[["adt.pca"]] <- do.call(runPca, c(list(store[["adt.normalized"]], block=block, num.threads=num.threads), runPca.args))
    } else {
        store["adt.pca"] <- list(NULL)
    }

    if (!is.null(crispr.x)) {
        store[["crispr.pca"]] <- do.call(runPca, c(list(store[["crispr.normalized"]], block=block, num.threads=num.threads), runPca.args))
    } else {
        store["crispr.pca"] <- list(NULL)
    }

    ############ Combining modalities and batches („• ᴗ •„) #############

    embeddings <- character()
    if (use.rna.pcs && !is.null(rna.x)) {
        embeddings <- c(embeddings, "rna.pca")
    }
    if (use.adt.pcs && !is.null(adt.x)) {
        embeddings <- c(embeddings, "adt.pca")
    }
    if (use.crispr.pcs && !is.null(crispr.x)) {
        embeddings <- c(embeddings, "crispr.pca")
    }

    if (length(embeddings) == 0) {
        stop("at least one 'use.*' must be True")
    }
    if (length(embeddings) == 1) {
        store[["combined.pca"]] <- embeddings[1]
        chosen.pcs <- store[[embeddings[1]]]$components
    } else {
        embeddings <- lapply(embeddings, function(e) store[[e]]$components)
        store[["combined.pca"]] <- do.call(scaleByNeighbors, c(list(embeddings, BNPARAM=BNPARAM, num.threads=num.threads), scaleByNeighbors.args))
        chosen.pcs <- store[["combined.pca"]]$combined
    }

    store[["block"]] <- block
    if (is.null(block)) {
        store["mnn.corrected"] <- list(NULL)
    } else {
        store[["mnn.corrected"]] <- do.call(correctMnn, c(list(chosen.pcs, block=block, BNPARAM=BNPARAM, num.threads=num.threads), correctMnn.args))
        chosen.pcs <- store[["mnn.corrected"]]$corrected
    }

    ############ Assorted neighbor-related stuff ⸜(⸝⸝⸝´꒳`⸝⸝⸝)⸝ #############

    nn.out <- do.call(runAllNeighborSteps, c(
        list(
            chosen.pcs,
            runUmap.args=runUmap.args, 
            runTsne.args=runTsne.args, 
            buildSnnGraph.args=buildSnnGraph.args,
            clusterGraph.args=clusterGraph.args,
            BNPARAM=BNPARAM,
            num.threads=num.threads
        ),
        runAllNeighborSteps.args
    ))

    store["tsne"] <- nn.out["runTsne"]
    store["umap"] <- nn.out["runUmap"]
    store["snn.graph"] <- nn.out["buildSnnGraph"]
    store["graph.clusters"] <- nn.out["clusterGraph"]

    ############ Finding markers (˶˃ ᵕ ˂˶) #############

    if (!is.null(kmeans.clusters)) {
        store[["kmeans.clusters"]] <- do.call(clusterKmeans, c(list(chosen.pcs, k=kmeans.clusters, num.threads=num.threads), clusterKmeans.args))
    } else {
        store["kmeans.clusters"] <- list(NULL)
    }

    store["clusters"] <- list(NULL)
    for (c in clusters.for.markers) {
        if (c == "graph") {
            candidate <- store[["graph.clusters"]]
            if (!is.null(candidate)) {
                store[["clusters"]] <- candidate$membership
                break
            }
        } else if (c == "kmeans") {
            candidate <- store[["kmeans.clusters"]]
            if (!is.null(candidate)) {
                store[["clusters"]] <- candidate$clusters
                break
            }
        }
    }

    store["rna.markers"] <- list(NULL)
    store["adt.markers"] <- list(NULL)
    store["crispr.markers"] <- list(NULL)

    clusters <- store[["clusters"]]
    if (!is.null(clusters)) {
        if (!is.null(rna.x)) {
            store[["rna.markers"]] <- do.call(scoreMarkers, c(list(store[["rna.normalized"]], groups=clusters, num.threads=num.threads, block=block), scoreMarkers.args))
        }
        if (!is.null(adt.x)) {
            store[["adt.markers"]] <- do.call(scoreMarkers, c(list(store[["adt.normalized"]], groups=clusters, num.threads=num.threads, block=block), scoreMarkers.args))
        }
        if (!is.null(crispr.x)) {
            store[["crispr.markers"]] <- do.call(scoreMarkers, c(list(store[["crispr.normalized"]], groups=clusters, num.threads=num.threads, block=block), scoreMarkers.args))
        }
    }

    store
}

convertQcMetricsToDataFrame <- function(metrics, flatten = TRUE) {
    nrow <- length(metrics$sum) # all QC functions produce 'sum'.
    output <- S4Vectors::make_zero_col_DFrame(nrow=nrow)

    for (n in names(metrics)) {
        current <- metrics[[n]]
        if (!is.list(current)) {
            output[[n]] <- current
        } else {
            if (is.null(names(current))) {
                names(current) <- seq_along(current)
            }
            if (flatten) {
                for (x in names(current)) {
                    output[[paste0(n, ".", x)]] <- current[[x]]
                }
            } else {
                tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow)
                for (x in names(current)) {
                    tmp[[x]] <- current[[x]]
                }
                output[[n]] <- tmp
            }
        }
    }

    output
}

#' Convert analysis results into a SingleCellExperiment
#'
#' Convert results from \code{\link{analyze}} into a \link[SingleCellExperiment]{SingleCellExperiment} for further analysis with Bioconductor packages.
#'
#' @param results List of results produced by \code{\link{analyze}}.
#' @param main.modality String specifying the modality to use as the main experiment of a \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param flatten.qc.subsets Logical scalar indicating whether QC metrics for subsets should be flattened in the column data.
#' If \code{FALSE}, subset metrics are reported as a nested \link[S4Vectors]{DataFrame}.
#' @param include.per.block.variances Logical scalar indicating whether the per-block variances should be reported as a nested \link[S4Vectors]{DataFrame} in the row data.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing most of the analysis results.
#' Filtered and normalized matrices are stored in the assays.
#' QC metrics, size factors and clusterings are stored in the column data.
#' Gene variances are stored in the row data.
#' PCA, t-SNE and UMAP results are stored in the reduced dimensions.
#' Further modalities are stored as alternative experiments.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{analyze}}, to generate \code{results}.
#'
#' @export
convertAnalyzeResults <- function(
    results,
    main.modality = NULL,
    flatten.qc.subsets = TRUE,
    include.per.block.variances = FALSE)
{
    objects <- list()

    maybe.filter <- function(x, ncells) {
        if (nrow(x) > ncells) {
            x[results$combined.qc.filter,,drop=FALSE]
        } else {
            x
        }
    }

    if (!is.null(results$rna.filtered)) {
        rna.sce <- SingleCellExperiment::SingleCellExperiment(
            assays=list(filtered=results$rna.filtered, normalized=results$rna.normalized),
            reducedDims=list(pca=t(results$rna.pca$components))
        )

        cd <- SummarizedExperiment::colData(rna.sce)
        qcdf <- convertQcMetricsToDataFrame(results$rna.qc.metrics, flatten=flatten.qc.subsets)
        cd <- cbind(cd, maybe.filter(qcdf, nrow(cd)))
        SummarizedExperiment::colData(rna.sce) <- cd
        SingleCellExperiment::sizeFactors(rna.sce) <- results$rna.size.factors

        rd <- SummarizedExperiment::rowData(rna.sce)
        rd <- cbind(rd, results$rna.gene.variances$statistics)
        if (include.per.block.variances) {
            per.block <- results$rna.gene.variances$per.block
            if (!is.null(per.block)) {
                for (i in seq_along(per.block)) {
                    per.block[[i]] <- S4Vectors::I(S4Vectors::DataFrame(per.block[[i]]))
                }
                rd$per.block <- S4Vectors::DataFrame(per.block)
            }
        }
        hvgs <- logical(nrow(rd))
        hvgs[results$rna.highly.variable.genes] <- TRUE
        rd$is.highly.variable <- hvgs
        SummarizedExperiment::rowData(rna.sce) <- rd

        objects$rna <- rna.sce
    }

    if (!is.null(results$adt.filtered)) {
        adt.sce <- SingleCellExperiment::SingleCellExperiment(
            assays=list(filtered=results$adt.filtered, normalized=results$adt.normalized),
            reducedDims=list(pca=t(results$adt.pca$components))
        )

        cd <- SummarizedExperiment::colData(adt.sce)
        qcdf <- convertQcMetricsToDataFrame(results$adt.qc.metrics, flatten=flatten.qc.subsets)
        cd <- cbind(cd, maybe.filter(qcdf, nrow(cd)))
        SummarizedExperiment::colData(adt.sce) <- cd
        SingleCellExperiment::sizeFactors(adt.sce) <- results$adt.size.factors

        objects$adt <- adt.sce
    }

    if (!is.null(results$crispr.filtered)) {
        crispr.sce <- SingleCellExperiment::SingleCellExperiment(
            assays=list(filtered=results$crispr.filtered, normalized=results$crispr.normalized),
            reducedDims=list(pca=t(results$crispr.pca$components))
        )

        cd <- SummarizedExperiment::colData(crispr.sce)
        qcdf <- convertQcMetricsToDataFrame(results$crispr.qc.metrics)
        cd <- cbind(cd, maybe.filter(qcdf, nrow(cd)))
        SummarizedExperiment::colData(crispr.sce) <- cd
        SingleCellExperiment::sizeFactors(crispr.sce) <- results$crispr.size.factors

        objects$crispr <- crispr.sce
    }

    if (is.null(main.modality)) {
        main.modality <- names(objects)[1]
    }

    sce <- objects[[main.modality]]
    SingleCellExperiment::mainExpName(sce) <- main.modality
    objects[[main.modality]] <- NULL
    SingleCellExperiment::altExps(sce) <- objects

    if (!is.character(results$combined.pca)) {
        SingleCellExperiment::reducedDim(sce, "combined.pca") <- t(results$combined.pca$combined)
    }

    if (!is.null(results$mnn.corrected)) {
        sce$block <- results$block
        SingleCellExperiment::reducedDim(sce, "mnn.corrected") <- t(results$mnn.corrected$corrected)
    }

    if (!is.null(results$tsne)) {
        SingleCellExperiment::reducedDim(sce, "tsne") <- results$tsne
    }

    if (!is.null(results$umap)) {
        SingleCellExperiment::reducedDim(sce, "umap") <- results$umap
    }

    if (!is.null(results$clusters)) {
        sce$clusters <- results$clusters
    }

    sce
}
