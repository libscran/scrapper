#' Scale and combine multiple embeddings in a SingleCellExperiment
#'
#' Scale embeddings for different modalities to equalize their intra-population variance, and combine them into a single embedding for downstream analysis.
#' This calls \code{\link{scaleByNeighbors}} on the reduced dimensions of the main/alternative experiments in a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param altexp.reddims Named list of character or integer vectors.
#' Each entry is named after an alternative experiment.
#' Each vector contains the names/indices of the \code{\link[SingleCellExperiment]{reducedDim}} embeddings from that experiment to be combined.
#' @param main.reddims Character or integer vector specifying the names/indices of the \code{\link[SingleCellExperiment]{reducedDim}} entries from \code{x} to be combined.
#' @param num.neighbors Number of neighbors used to define the scaling factor, passed to \code{\link{scaleByNeighbors}}.
#' @param block Block assignment for each cell, passed to \code{\link{scaleByNeighbors}}.
#' @param num.threads Number of threads for the neighbor search, passed to \code{\link{scaleByNeighbors}}.
#' @param BNPARAM Algorithm for the nearest neighbor search, passed to \code{\link{scaleByNeighbors}}.
#' @param more.scale.args Named list of additional arguments to pass to \code{\link[scrapper]{scaleByNeighbors}}.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the combined embeddings.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store additional metrics.
#' If \code{NULL}, additional metrics are not stored.
#' @param delayed.transpose Logical scalar indicating whether to delay the transposition when storing coordinates in the \code{\link[SingleCellExperiment]{reducedDims}}.
#'
#' @return \code{x} is returned with the combined embeddings stored in its \code{rowData}.
#' The scaling factors for all embeddings are stored in the \code{metadata}.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestAdtData.se("pca")
#' sce <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT="PCA"))
#' reducedDimNames(sce) 
#' metadata(sce)$combined
#'
#' @export
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom S4Vectors metadata metadata<-
scaleByNeighbors.se <- function(
    x,
    altexp.reddims,
    main.reddims = "PCA", 
    num.neighbors = 20,
    block = NULL,
    BNPARAM = AnnoyParam(),
    num.threads = 1,
    more.scale.args = list(),
    output.name = "combined",
    meta.name = "combined",
    delayed.transpose = FALSE
) {
    all.embeddings <- list()
    main.reddims <- .sanitizeReddims(main.reddims, SingleCellExperiment::reducedDimNames(x))
    for (r in main.reddims) {
        all.embeddings <- append(all.embeddings, list(.getTransposedReddim(x, r)))
    }

    altexp.reddims <- altexp.reddims[!duplicated(names(altexp.reddims))]
    for (ae in names(altexp.reddims)) {
        ae.se <- SingleCellExperiment::altExp(x, ae)
        cur.reddim <- .sanitizeReddims(altexp.reddims[[ae]], SingleCellExperiment::reducedDimNames(ae.se))
        for (r in cur.reddim) {
            all.embeddings <- append(all.embeddings, list(.getTransposedReddim(ae.se, r)))
        }
        altexp.reddims[[ae]] <- cur.reddim
    }

    out <- .call(
        scaleByNeighbors,
        list(all.embeddings),
        list(num.neighbors=num.neighbors, block=block, num.threads=num.threads, BNPARAM=BNPARAM),
        more.scale.args
    )

    x <- .addTransposedReddim(x, output.name, out$combined, delayed.transpose)
    if (!is.null(meta.name)) {
        out$combined <- NULL

        # Formatting it in the same manner as the arguments.
        counter <- 1L
        out$main.scaling <- numeric(0)
        for (r in main.reddims) {
            out$main.scaling[[r]] <- out$scaling[[counter]]
            counter <- counter + 1L
        }

        out$altexp.scaling <- list()
        for (ae in names(altexp.reddims)) {
            current <- numeric(0)
            for (r in altexp.reddims[[ae]]) {
                current[[r]] <- out$scaling[[counter]]
                counter <- counter + 1L
            }
            out$altexp.scaling[[ae]] <- current
        }

        out$scaling <- NULL
        metadata(x)[[meta.name]] <- out
    }

    x
}

.sanitizeReddims <- function(reddims, all.reddims) {
    reddims <- unique(reddims)
    if (is.integer(reddims)) {
        all.reddims[reddims]
    } else {
        reddims
    }
}
