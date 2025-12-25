#' t-SNE on a SummarizedExperiment
#'
#' Generate a t-SNE visualization from an existing embedding,
#' by calling \code{\link{runUmap}} on a reduced dimension entry in \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param perplexity Perplexity to use in the t-SNE algorithm, passed to \code{\link{runTsne}}.
#' @param num.threads Number of threads for the neighbor search and optimization, passed to \code{\link{runTsne}}.
#' @param more.tsne.args Named list of further arguments to pass to \code{\link{runTsne}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#' @param output.name String containing the name of the output \code{\link[SingleCellExperiment]{reducedDim}}. 
#'
#' @return \code{x} is returned with the t-SNE coordinates stored in the \code{reducedDim}.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("pca")
#' # Using fewer iterations for a faster-running example.
#' sce <- runTsne.se(sce, more.tsne.args=list(max.iterations=50))
#' head(reducedDim(sce, "TSNE"))
#'
#' @export
runTsne.se <- function(
    x,
    perplexity = 30,
    num.threads = 1,
    more.tsne.args = list(),
    reddim.type = "PCA", 
    output.name = "TSNE"
) {
    out <- .call(
        runTsne,
        list(.getTransposedReddim(x, reddim.type)),
        list(perplexity=perplexity, num.threads=num.threads),
        more.tsne.args
    )

    .addTsneResults(x, output.name, out)
}

.addTsneResults <- function(x, output.name, res) {
    SingleCellExperiment::reducedDim(x, output.name) <- res
    x
}
