#' MNN correction on a SingleCellExperiment
#'
#' Correct batch effects from an existing embedding with mutual nearest neighbors (MNNs),
#' by calling \code{\link{correctMnn}} on a reduced dimension entry of a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param block Block assignment for each cell, passed to \code{\link{correctMnn}}.
#' @param BNPARAM Algorithm for the nearest neighbor search, passed to \code{\link{correctMnn}}.
#' @param num.threads Number of threads, passed to \code{\link{correctMnn}}.
#' @param more.mnn.args Named list of additional arguments to pass to \code{\link{correctMnn}}.
#' @param reddim.type String or integer specifying the \code{\link[SingleCellExperiment]{reducedDim}} entry on which to perform MNN correction.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the corrected embedding.
#' @param delayed.transpose Logical scalar indicating whether to delay the transposition when storing coordinates in the \code{\link[SingleCellExperiment]{reducedDims}}.
#'
#' @return \code{x} is returned with the corrected embedding stored as a \code{reducedDim} entry.
#' @author Aaron Lun
#'
#' @examples
#' sce <- getTestRnaData.se("pca")
#' # Treating the tissue of origin as the batch.
#' sce <- correctMnn.se(sce, sce$tissue)
#' reducedDimNames(sce)
#'
#' @export
#' @importFrom BiocNeighbors AnnoyParam
correctMnn.se <- function(
    x,
    block,
    BNPARAM = AnnoyParam(),
    num.threads = 1,
    more.mnn.args = list(),
    reddim.type = "PCA", 
    output.name = "MNN",
    delayed.transpose = FALSE
) {
    out <- .call(
        correctMnn,
        list(.getTransposedReddim(x, reddim.type)),
        list(block=block, num.threads=num.threads, BNPARAM=BNPARAM),
        more.mnn.args
    )

    .addTransposedReddim(x, output.name, out$corrected, delayed.transpose)
}
