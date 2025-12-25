#' UMAP on a SummarizedExperiment
#'
#' Generate a UMAP visualization from an existing embedding,
#' by calling \code{\link{runUmap}} on a reduced dimension entry in \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param num.dim Number of dimensions in the output embedding, passed to \code{\link{runUmap}}.
#' @param num.neighbors Number of neighbors for constructing the fuzzy sets, passed to \code{\link{runUmap}}.
#' @param min.dist Minimum distance between observations, passed to \code{\link{runUmap}}.
#' @param num.threads Number of threads for the UMAP, passed to \code{\link{runUmap}}.
#' @param more.umap.args Named list of further arguments to pass to \code{\link[scrapper]{runUmap}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#' @param output.name String containing the name of the output \code{\link[SingleCellExperiment]{reducedDim}}. 
#'
#' @return \code{x} is returned with the UMAP coordinates stored in the \code{reducedDim}.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("pca")
#' # Using fewer epochs for a faster-running example.
#' sce <- runUmap.se(sce, more.umap.args=list(num.epochs=50))
#' head(reducedDim(sce, "UMAP"))
#'  
#' @export
runUmap.se <- function(
    x,
    num.dim = 2,
    min.dist = 0.1,
    num.neighbors = 15,
    num.threads = 1,
    more.umap.args = list(),
    reddim.type = "PCA", 
    output.name = "UMAP"
) {
    res <- .call(
        runUmap,
        list(.getTransposedReddim(x, reddim.type)),
        list(num.dim=num.dim, min.dist=min.dist, num.neighbors=num.neighbors, num.threads=num.threads),
        more.umap.args
    )

    .addUmapResults(x, output.name, res)
}

.addUmapResults <- function(x, output.name, res) {
    SingleCellExperiment::reducedDim(x, output.name) <- res
    x
}
