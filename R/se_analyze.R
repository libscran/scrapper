#' Analyze single-cell data from a SummarizedExperiment
#'
#' Execute a simple single-cell analysis pipeline, starting from a count matrix and ending with clusters, visualizations and markers.
#' This also supports integration of multiple modalities and correction of batch effects.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param rna.altexp String or integer specifying the alternative experiment of \code{x} containing the RNA data.
#' If \code{NA}, the main experiment is assumed to contain the RNA data.
#' If \code{NULL}, it is assumed that no RNA data is available.
#' @param adt.altexp String or integer specifying the alternative experiment of \code{x} containing the ADT data.
#' If \code{NA}, the main experiment is assumed to contain the ADT data.
#' If \code{NULL}, it is assumed that no ADT data is available.
#' @param crispr.altexp String or integer specifying the alternative experiment of \code{x} containing the CRISPR data.
#' If \code{NA}, the main experiment is assumed to contain the CRISPR data.
#' If \code{NULL}, it is assumed that no CRISPR data is available.
#' @param rna.assay.type String or integer specifying the assay containing the RNA count data.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param adt.assay.type String or integer specifying the assay containing the ADT count data.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param crispr.assay.type String or integer specifying the assay containing the CRISPR count data.
#' Only used if \code{crispr.altexp} is not \code{NULL}.
#' @param block Vector or factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' Alternatively \code{NULL}, if all cells are from the same block.
#' @param block.name String specifying the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the blocking factor.
#' Only used if \code{block} is not \code{NULL}.
#' If \code{NULL}, the blocking factor is not stored in the \code{colData}.
#' @param rna.qc.subsets Passed to \code{\link{quickRnaQc.se}} as the \code{subsets} argument.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param rna.qc.output.prefix Passed to \code{\link{quickRnaQc.se}} as the \code{output.prefix} argument.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param more.rna.qc.args Named list of additional arguments to pass to \code{\link{quickRnaQc.se}}.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param adt.qc.subsets Passed to \code{\link{quickAdtQc.se}} as the \code{subsets} argument.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param adt.qc.output.prefix Passed to \code{\link{quickAdtQc.se}} as the \code{output.prefix} argument.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param more.adt.qc.args Named list of additional arguments to pass to \code{\link{quickAdtQc.se}}.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param crispr.qc.output.prefix Passed to \code{\link{quickCrisprQc.se}} as the \code{output.prefix} argument.
#' Only used if \code{crispr.altexp} is not \code{NULL}.
#' @param more.crispr.qc.args Named list of additional arguments to pass to \code{\link{quickCrisprQc.se}}.
#' Only used if \code{crispr.altexp} is not \code{NULL}.
#' @param filter.cells Logical scalar indicating whether to filter \code{x} to only retain high-quality cells in all modalities.
#' If \code{FALSE}, QC metrics and thresholds are still computed but are not used to filter the count matrices.
#' @param rna.norm.output.name Passed to \code{\link{normalizeRnaCounts.se}} as the \code{output.name} argument.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param more.rna.norm.args Named list of arguments to pass to \code{\link{normalizeRnaCounts.se}}.
#' Only used if \code{rna.altexp} is not \code{NULL}.
#' @param adt.norm.output.name Passed to \code{\link{normalizeAdtCounts.se}} as the \code{output.name} argument.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param more.adt.norm.args Named list of arguments to pass to \code{\link{normalizeAdtCounts.se}}.
#' Only used if \code{adt.altexp} is not \code{NULL}.
#' @param crispr.norm.output.name Passed to \code{\link{normalizeCrisprCounts.se}} as the \code{output.name} argument.
#' Only used if \code{crispr.altexp} is not \code{NULL}.
#' @param more.crispr.norm.args Named list of arguments to pass to \code{\link{normalizeCrisprCounts.se}}.
#' Only used if \code{crispr.altexp} is not \code{NULL}.
#' @param rna.hvg.output.prefix Passed to \code{\link{chooseRnaHvgs.se}} as the \code{output.prefix} argument.
#' Only used if \code{rna.altexp} is provided.
#' @param more.rna.hvg.args Named list of arguments to pass to \code{\link{chooseRnaHvgs.se}}.
#' Only used if \code{rna.altexp} is provided.
#' @param rna.pca.output.name Passed to \code{\link{runPca.se}} as the \code{output.name} argument.
#' Only used if \code{rna.altexp} is provided.
#' @param more.rna.pca.args Named list of arguments to pass to \code{\link{runPca.se}}.
#' Only used if \code{rna.altexp} is provided.
#' @param adt.pca.output.name Passed to \code{\link{runPca.se}} as the \code{output.name} argument.
#' Only used if \code{adt.altexp} is provided.
#' @param more.adt.pca.args Named list of arguments to pass to \code{\link{runPca.se}}.
#' Only used if \code{adt.altexp} is provided.
#' @param use.rna.pcs Logical scalar indicating whether to use the RNA-derived PCs for downstream steps (i.e., clustering, visualization).
#' Only used if \code{rna.altexp} is provided.
#' @param use.adt.pcs Logical scalar indicating whether to use the ADT-derived PCs for downstream steps (i.e., clustering, visualization).
#' Only used if \code{adt.altexp} is provided.
#' @param scale.output.name Passed to \code{\link{scaleByNeighbors.se}} as the \code{output.name} argument.
#' Only used if multiple modalities are available and their corresponding \code{use.*.pcs} arguments are \code{TRUE}.
#' @param more.scale.args Named list of arguments to pass to \code{\link{scaleByNeighbors.se}}.
#' Only used if multiple modalities are available and their corresponding \code{use.*.pcs} arguments are \code{TRUE}.
#' @param mnn.output.name Passed to \code{\link{correctMnn.se}} as the \code{output.name} argument.
#' Only used if \code{block} is supplied.
#' @param more.mnn.args Named list of arguments to pass to \code{\link{correctMnn.se}}.
#' Only used if \code{block} is supplied.
#' @param more.tsne.args Passed to \code{\link{runAllNeighborSteps.se}}.
#' @param more.umap.args Passed to \code{\link{runAllNeighborSteps.se}}.
#' @param more.build.graph.args Passed to \code{\link{runAllNeighborSteps.se}}. 
#' @param cluster.graph.output.name Passed to \code{\link{runAllNeighborSteps.se}} as \code{cluster.output.name}.
#' @param more.cluster.graph.args Passed to \code{\link{runAllNeighborSteps.se}}. 
#' @param more.neighbor.args Named list of arguments to pass to \code{\link{runAllNeighborSteps.se}}.
#' @param kmeans.clusters Passed to \code{\link{clusterKmeans.se}} as the \code{k} argument.
#' If \code{NULL}, k-means clustering is not performed.
#' @param kmeans.clusters.output.name Passed to \code{\link{clusterKmeans.se}} as the \code{output.name} argument.
#' Ignored if \code{kmeans.clusters = NULL}.
#' @param more.kmeans.args Named list of arguments to pass to \code{\link{clusterKmeans.se}}.
#' Ignored if \code{kmeans.clusters = NULL}.
#' @param clusters.for.markers Character vector of clustering algorithms (either \code{"graph"} or \code{"kmeans"}, specifying the clustering to be used for marker detection.
#' The first available clustering will be chosen.
#' @param more.rna.marker.args Named list of arguments to pass to \code{\link{scoreMarkers.se}} for the RNA data.
#' Ignored if no suitable clusterings are available or if \code{rna.altexp=NULL}.
#' @param more.adt.marker.args Named list of arguments to pass to \code{\link{scoreMarkers.se}} for the ADT data.
#' Ignored if no suitable clusterings are available or if \code{adt.altexp=NULL}.
#' @param more.crispr.marker.args Named list of arguments to pass to \code{\link{scoreMarkers.se}} for the CRISPR data.
#' Ignored if no suitable clusterings are available or if \code{crispr.altexp=NULL}.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} instance specifying the nearest-neighbor search algorithm to use.
#' @param num.threads Integer scalar specifying the number of threads to use in each step.
#'
#' @return \code{x} is returned with the results of the analysis.
#'
#' @details
#' This function is equivalent to:
#' \itemize{
#' \item Running \code{\link{quickRnaQc.se}}, \code{\link{quickAdtQc.se}} and/or \code{\link{quickCrisprQc.se}}, for quality control.
#' \item Subsetting \code{x} to only retain the high-quality cells in all modalities, based on \code{filter.cells}.
#' \item Running \code{\link{normalizeRnaCounts.se}}, \code{\link{normalizeAdtCounts.se}} and/or \code{\link{normalizeCrisprCounts.se}}, for normalization.
#' \item Running \code{\link{chooseRnaHvgs.se}} to identify highly variable genes.
#' \item Running \code{\link{runPca.se}} on the RNA and/or ADT data.
#' \item Running \code{\link{scaleByNeighbors.se}} if multiple modalities are present.
#' \item Running \code{\link{correctMnn.se}} if multiple batches are present in \code{block}.
#' \item Running \code{\link{runAllNeighborSteps.se}} to obtain t-SNE and UMAP coordinates, and to perform graph-based clustering.
#' \item Running \code{\link{clusterKmeans.se}} to perform k-means clustering, depending on \code{kmeans.clusters}. 
#' \item Running \code{\link{scoreMarkers.se}} to compute markers for each modality based on one of the clusterings. 
#' }
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("start")
#' res <- analyze.se(
#'     sce, 
#'     rna.qc.subsets=list(mito=grep("^mt-", rownames(sce))),
#'     num.threads=2 # keep R CMD check happy
#' )
#' assayNames(res$x)
#' reducedDimNames(res$x)
#' colData(res$x)
#' previewMarkers(res$markers$rna[[1]], "cohens.d.mean")
#'
#' @export
#' @importFrom methods is
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom DelayedArray DelayedArray
analyze.se <- function(
    x,
    rna.altexp = NA,
    adt.altexp = NULL,
    crispr.altexp = NULL,
    rna.assay.type = "counts",
    adt.assay.type = "counts",
    crispr.assay.type = "counts",
    block = NULL,
    block.name = "block",
    rna.qc.subsets = list(),
    rna.qc.output.prefix = NULL,
    more.rna.qc.args = list(),
    adt.qc.output.prefix = NULL,
    adt.qc.subsets = list(),
    more.adt.qc.args = list(),
    crispr.qc.output.prefix = NULL,
    more.crispr.qc.args = list(),
    filter.cells = TRUE,
    rna.norm.output.name = "logcounts",
    more.rna.norm.args = list(),
    adt.norm.output.name = "logcounts",
    more.adt.norm.args = list(),
    crispr.norm.output.name = "logcounts",
    more.crispr.norm.args = list(),
    rna.hvg.output.prefix = NULL,
    more.rna.hvg.args = list(),
    rna.pca.output.name = "PCA",
    more.rna.pca.args = list(),
    adt.pca.output.name = "PCA",
    more.adt.pca.args = list(),
    use.rna.pcs = TRUE,
    use.adt.pcs = TRUE,
    scale.output.name = "combined",
    more.scale.args = list(),
    mnn.output.name = "MNN",
    more.mnn.args = list(),
    more.umap.args = list(),
    more.tsne.args = list(),
    cluster.graph.output.name = "graph.cluster",
    more.build.graph.args = list(),
    more.cluster.graph.args = list(),
    more.neighbor.args = list(),
    kmeans.clusters = NULL,
    kmeans.clusters.output.name = "kmeans.cluster",
    more.kmeans.args = list(),
    clusters.for.markers = c("graph", "kmeans"),
    more.rna.marker.args = list(),
    more.adt.marker.args = list(),
    more.crispr.marker.args = list(),
    BNPARAM = AnnoyParam(),
    num.threads = 3L
) {
    store <- list()

    ############ Quality control o(*°▽°*)o #############

    collected.filters <- list()

    if (!is.null(rna.altexp)) {
        tmp <- .getModality(x, rna.altexp)
        tmp <- .call(
            quickRnaQc.se, 
            list(tmp, output.prefix=rna.qc.output.prefix),
            list(subsets=rna.qc.subsets, block=block, num.threads=num.threads, assay.type=rna.assay.type),
            more.rna.qc.args
        )
        collected.filters$rna <- .extractOrError(SummarizedExperiment::colData(tmp), paste0(rna.qc.output.prefix, "keep"))
        x <- .assignModality(x, rna.altexp, tmp)
    }

    if (!is.null(adt.altexp)) {
        tmp <- .getModality(x, adt.altexp)
        tmp <- .call(
            quickAdtQc.se,
            list(tmp, output.prefix=adt.qc.output.prefix),
            list(subsets=adt.qc.subsets, block=block, num.threads=num.threads, assay.type=adt.assay.type),
            more.adt.qc.args
        )
        collected.filters$adt <- .extractOrError(SummarizedExperiment::colData(tmp), paste0(adt.qc.output.prefix, "keep"))
        x <- .assignModality(x, adt.altexp, tmp)
    }

    if (!is.null(crispr.altexp)) {
        tmp <- .getModality(x, crispr.altexp)
        tmp <- .call(
            quickCrisprQc.se, 
            list(tmp, output.prefix=crispr.qc.output.prefix),
            list(num.threads=num.threads, block=block, assay.type=crispr.assay.type),
            more.crispr.qc.args
        )
        collected.filters$crispr <- .extractOrError(SummarizedExperiment::colData(tmp), paste0(crispr.qc.output.prefix, "keep"))
        x <- .assignModality(x, crispr.altexp, tmp)
    }

    # Combining all filters.
    combined.qc.filter <- !logical(ncol(x))
    for (mod in c("rna", "adt", "crispr")) {
        keep <- collected.filters[[mod]]
        if (!is.null(keep)) {
            combined.qc.filter <- combined.qc.filter & keep
        }
    }
    if (length(collected.filters) > 1) {
        SummarizedExperiment::colData(x)[["combined.keep"]] <- combined.qc.filter
    }

    if (filter.cells) {
        x <- .delayifyAssays(x)
        x <- x[,combined.qc.filter]
        if (!is.null(block)) {
            block <- block[combined.qc.filter]
        }
    }

    ############ Normalization ( ꈍᴗꈍ) #############

    if (!is.null(rna.altexp)) {
        tmp <- .getModality(x, rna.altexp)
        lib.sizes <- .extractOrError(SummarizedExperiment::colData(tmp), paste0(rna.qc.output.prefix, "sum"))
        tmp <- .call(
            normalizeRnaCounts.se,
            list(tmp, assay.type=rna.assay.type, output.name=rna.norm.output.name),
            list(size.factors=lib.sizes, block=block),
            more.rna.norm.args
        )
        x <- .assignModality(x, rna.altexp, tmp)
    }

    if (!is.null(adt.altexp)) {
        tmp <- .getModality(x, adt.altexp)
        tmp <- .call(
            normalizeAdtCounts.se, 
            list(tmp, assay.type=adt.assay.type, output.name=adt.norm.output.name),
            list(block=block),
            more.adt.norm.args
        )
        x <- .assignModality(x, adt.altexp, tmp)
    }

    if (!is.null(crispr.altexp)) {
        tmp <- .getModality(x, crispr.altexp)
        lib.sizes <- .extractOrError(SummarizedExperiment::colData(tmp), paste0(crispr.qc.output.prefix, "sum"))
        tmp <- .call(
            normalizeCrisprCounts.se,
            list(tmp, assay.type=crispr.assay.type, output.name=crispr.norm.output.name),
            list(size.factors=lib.sizes, block=block),
            more.crispr.norm.args
        )
        x <- .assignModality(x, crispr.altexp, tmp)
    }

    ############ Variance modelling (～￣▽￣)～ #############

    if (!is.null(rna.altexp)) {
        tmp <- .getModality(x, rna.altexp)
        tmp <- .call(
            chooseRnaHvgs.se,
            list(tmp, output.prefix=rna.hvg.output.prefix),
            list(block=block),
            more.rna.hvg.args
        )
        x <- .assignModality(x, rna.altexp, tmp)
    }

    ############ Principal components analysis \(>⩊<)/ #############

    if (!is.null(rna.altexp)) {
        tmp <- .getModality(x, rna.altexp)
        is.hvg <- .extractOrError(SummarizedExperiment::rowData(tmp), paste0(rna.hvg.output.prefix, "hvg"))
        tmp <- .call(
            runPca.se,
            list(tmp, assay.type=rna.norm.output.name, output.name=rna.pca.output.name),
            list(features=is.hvg, block=block, delayed.transpose=TRUE),
            more.rna.pca.args
        )
        x <- .assignModality(x, rna.altexp, tmp)
    }

    if (!is.null(adt.altexp)) {
        tmp <- .getModality(x, adt.altexp)
        tmp <- .call(
            runPca.se,
            list(tmp, assay.type=adt.norm.output.name, output.name=adt.pca.output.name),
            list(features=NULL, block=block, delayed.transpose=TRUE),
            more.adt.pca.args
        )
        x <- .assignModality(x, adt.altexp, tmp)
    }

    ############ Combining modalities („• ᴗ •„) #############

    embeddings <- character()
    if (use.rna.pcs && !is.null(rna.altexp)) {
        embeddings <- c(embeddings, "rna")
    }
    if (use.adt.pcs && !is.null(adt.altexp)) {
        embeddings <- c(embeddings, "adt")
    }

    target.embedding <- NULL
    if (length(embeddings) == 0) {
        stop("at least one 'use.*' must be true")

    } else if (length(embeddings) == 1L) {
        if (embeddings == "rna") {
            target.embedding <- .defineSingleTargetEmbedding(x, rna.altexp, rna.pca.output.name)
        } else if (embeddings == "adt") {
            target.embedding <- .defineSingleTargetEmbedding(x, adt.altexp, adt.pca.output.name)
        }

    } else if (length(embeddings) > 1) {
        all.reddims <- list(main=character(0), altexp=list())
        if ("rna" %in% embeddings) {
            all.reddims <- .addSourceEmbeddingToScale(x, rna.altexp, rna.pca.output.name, all.reddims)
        }
        if ("adt" %in% embeddings) {
            all.reddims <- .addSourceEmbeddingToScale(x, adt.altexp, adt.pca.output.name, all.reddims)
        }

        x <- .call(
            scaleByNeighbors.se,
            list(x, output.name=scale.output.name, main.reddims=all.reddims$main, altexp.reddims=all.reddims$altexp, delayed.transpose=TRUE),
            list(block=block, num.threads=num.threads, BNPARAM=BNPARAM),
            more.scale.args
        )
        target.embedding <- scale.output.name
    }

    ############ Performing batch correction ⸜(｡˃ ᵕ ˂ )⸝ #############

    if (!is.null(block)) {
        x <- .call(
            correctMnn.se,
            list(x, output.name=mnn.output.name, reddim.type=target.embedding, delayed.transpose=TRUE),
            list(block=block, num.threads=num.threads, BNPARAM=BNPARAM),
            more.mnn.args
        )
        target.embedding <- mnn.output.name

        if (!is.null(block.name)) {
            SummarizedExperiment::colData(x)[[block.name]] <- block
        }
    }

    ############ Assorted neighbor-related stuff ⸜(⸝⸝⸝´꒳`⸝⸝⸝)⸝ #############

    x <- .call(
        runAllNeighborSteps.se,
        list(x, reddim.type=target.embedding, cluster.output.name=cluster.graph.output.name),
        list(
            more.umap.args=more.umap.args, 
            more.tsne.args=more.tsne.args, 
            more.build.graph.args=more.build.graph.args,
            more.cluster.graph.args=more.cluster.graph.args,
            BNPARAM=BNPARAM,
            num.threads=num.threads
        ),
        more.neighbor.args
    )

    ############ Maybe some k-means clustering ⸜(⸝⸝⸝´꒳`⸝⸝⸝)⸝ #############

    if (!is.null(kmeans.clusters)) {
        x <- .call(
            clusterKmeans.se,
            list(x, reddim.type=target.embedding, output.name=kmeans.clusters.output.name),
            list(k=kmeans.clusters, num.threads=num.threads),
            more.kmeans.args
        )
    }

    clusters.for.markers <- match.arg(clusters.for.markers, several.ok=TRUE)
    chosen.clusters <- NULL
    for (c in clusters.for.markers) {
        if (c == "graph") {
            name <- cluster.graph.output.name
        } else if (c == "kmeans") {
            name <- kmeans.clusters.output.name
        }
        if (!is.null(name) && name %in% colnames(SummarizedExperiment::colData(x))) {
            chosen.clusters <- SummarizedExperiment::colData(x)[[name]]
            break
        }
    }

    ############ Finding markers (˶˃ ᵕ ˂˶) #############

    markers <- list()

    if (!is.null(chosen.clusters)) {
        if (!is.null(rna.altexp)) {
            markers$rna <- .call(
                scoreMarkers.se,
                list(.getModality(x, rna.altexp), assay.type=rna.norm.output.name),
                list(groups=chosen.clusters, num.threads=num.threads, block=block),
                more.rna.marker.args
            )
        }

        if (!is.null(adt.altexp)) {
            markers$adt <- .call(
                scoreMarkers.se,
                list(.getModality(x, adt.altexp), assay.type=adt.norm.output.name),
                list(groups=chosen.clusters, num.threads=num.threads, block=block),
                more.adt.marker.args
            )
        }

        if (!is.null(crispr.altexp)) {
            markers$crispr <- .call(
                scoreMarkers.se,
                list(.getModality(x, crispr.altexp), assay.type=crispr.norm.output.name),
                list(groups=chosen.clusters, num.threads=num.threads, block=block),
                more.crispr.marker.args
            )
        }
    }

    list(x=x, markers=markers)
}

#' @importFrom methods is
.delayifyAssays <- function(x) {
    for (i in seq_along(SummarizedExperiment::assayNames(x))) {
        SummarizedExperiment::assay(x, i, withDimnames=FALSE) <- DelayedArray(SummarizedExperiment::assay(x, i, withDimnames=FALSE))
    }
    if (is(x, "SingleCellExperiment")) {
        for (i in seq_along(SingleCellExperiment::altExpNames(x))) {
            SingleCellExperiment::altExp(x, i, withDimnames=FALSE) <- .delayifyAssays(SingleCellExperiment::altExp(x, i, withDimnames=FALSE))
        }
    }
    x
}

.extractOrError <- function(df, name) {
    val <- df[[name]]
    if (is.null(val)) {
        stop("no DataFrame column named '", name, "'")
    }
    val
}

.getModality <- function(x, altexp) {
    if (is.na(altexp)) {
        x
    } else {
        SingleCellExperiment::altExp(x, altexp)
    }
}

.assignModality <- function(x, altexp, y) {
    if (is.na(altexp)) {
        y
    } else {
        SingleCellExperiment::altExp(x, altexp) <- y
        x
    }
}

.defineSingleTargetEmbedding <- function(x, altexp, output.name) {
    if (!is.na(altexp)) {
        if (is.numeric(altexp)) {
            nm <- SingleCellExperiment::altExpNames(x)[altexp]
        } else {
            nm <- altexp
        }
        names(output.name) <- nm
    }
    output.name
}

.addSourceEmbeddingToScale <- function(x, altexp, output.name, collected) {
    if (is.na(altexp)) {
        collected$main <- append(collected$main, output.name)
    } else {
        if (is.numeric(altexp)) {
            nm <- SingleCellExperiment::altExpNames(x)[altexp]
        } else {
            nm <- altexp
        }
        collected$altexp[[nm]] <- output.name
    }
    collected
}
