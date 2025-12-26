#' Run all nearest neighbor steps on a SummarizedExperiment
#'
#' Concurrently run all steps involving a nearest-neighbor search (t-SNE, UMAP and graph-based clustering) using the same nearest-neighbor index,
#' by calling \code{\link{runAllNeighborSteps}} on a reduced dimension entry of a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param umap.output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry to store the UMAP coordinates.
#' If \code{NULL}, the UMAP is not computed.
#' @param more.umap.args Named list of additional arguments to pass to \code{\link{runAllNeighborSteps}} as \code{runUmap.args}.
#' @param tsne.output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry to store the t-SNE coordinates.
#' If \code{NULL}, the t-SNE is not computed.
#' @param more.tsne.args Named list of additional arguments to pass to \code{\link{runAllNeighborSteps}} as \code{runTsne.args}.
#' @param build.graph.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store the nearest neighbor graph.
#' If \code{NULL}, the graph is not stored.
#' @param more.build.graph.args Named list of additional arguments to pass to \code{\link{runAllNeighborSteps}} as \code{buildSnnGraph.args}.
#' @param cluster.output.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the cluster assignments.
#' If \code{NULL}, graph-based clustering is not performed.
#' @param cluster.meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store additional clustering outputs.
#' If \code{NULL}, these additional outputs are not stored.
#' @param more.cluster.graph.args Named list of additional arguments to pass to \code{\link{runAllNeighborSteps}} as \code{clusterGraph.args}.
#' @param BNPARAM,num.threads Arguments to pass to \code{\link{runAllNeighborSteps}}.
#' @param more.neighbor.args Named list of additional arguments to pass to \code{\link{runAllNeighborSteps}}.
#' @param reddim.type String or integer specifying the \code{\link[SingleCellExperiment]{reducedDim}} entry on which to perform a nearest neighbor search.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#'
#' @return \code{x} is returned with additional coordinates stored in its \code{\link[SingleCellExperiment]{reducedDims}}
#' and clustering output in its \code{\link[SummarizedExperiment]{colData}}.
#' Additional information may also be stored in its \code{\link[S4Vectors]{metadata}}.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("pca")
#' sce <- runAllNeighborSteps.se(
#'    sce,
#'    more.tsne.args=list(max.iterations=50),
#'    more.umap.args=list(num.epochs=50),
#'    num.threads=2 # to keep R CMD check happy
#' )
#' reducedDimNames(sce)
#' table(sce$clusters)
#' 
#' @export 
#' @importFrom BiocNeighbors AnnoyParam
runAllNeighborSteps.se <- function(
    x,
    umap.output.name = "UMAP",
    more.umap.args = list(),
    tsne.output.name = "TSNE",
    more.tsne.args = list(),
    build.graph.name = NULL,
    more.build.graph.args = list(),
    cluster.output.name = "clusters",
    cluster.meta.name = NULL,
    more.cluster.graph.args = list(),
    BNPARAM = AnnoyParam(),
    num.threads = 3,
    more.neighbor.args = list(),
    reddim.type = "PCA"
) {
    if (is.null(umap.output.name)) {
        more.umap.args <- NULL 
    }
    if (is.null(tsne.output.name)) {
        more.tsne.args <- NULL 
    }
    if (is.null(cluster.output.name)) {
        more.cluster.graph.args <- NULL 
    }
    more.neighbor.args$return.graph <- !is.null(build.graph.name)

    outputs <- .call(
        runAllNeighborSteps,
        list(.getTransposedReddim(x, reddim.type)),
        list(
            runUmap.args=more.umap.args,
            runTsne.args=more.tsne.args,
            buildSnnGraph.args=more.build.graph.args,
            clusterGraph.args=more.cluster.graph.args,
            BNPARAM = BNPARAM,
            num.threads=num.threads
        ),
        more.neighbor.args
    )

    if (!is.null(outputs$runTsne)) {
        x <- .addTsneResults(x, outputs$runTsne, output.name=tsne.output.name)
    }
    if (!is.null(outputs$runUmap)) {
        x <- .addUmapResults(x, outputs$runUmap, output.name=umap.output.name)
    }
    if (!is.null(outputs$clusterGraph)) {
        x <- .addBuildGraphResults(x, outputs$buildSnnGraph, graph.name=build.graph.name)
        x <- .addClusterGraphResults(x, outputs$clusterGraph, output.name=cluster.output.name, meta.name=cluster.meta.name)
    }

    x
}
