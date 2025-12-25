#' Graph-based clustering of cells in a SingleCellExperiment 
#'
#' Construct a shared-nearest neighbor (SNN) graph from an existing low-dimensional embedding
#' by calling \code{\link{buildSnnGraph}} on a reduced dimension entry in a \link[SingleCellExperiment]{SingleCellExperiment}.
#' Then, apply community detection algorithms to obtain clusters of cells with \code{\link{clusterGraph}}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param num.neighbors Number of neighbors for constructing the graph, passed to \code{\link{buildSnnGraph}}. 
#' @param num.threads Number of threads for graph construction, passed to \code{\link{buildSnnGraph}}.
#' @param more.build.args Named list of further arguments to be passed to \code{\link{buildSnnGraph}}.
#' @param method Clustering method to use, passed to \code{\link{clusterGraph}}. 
#' @param resolution Resolution for the community detection method in \code{\link{clusterGraph}}.
#' This is either passed to \code{multilevel.resolution} or \code{leiden.resolution} depending on \code{method}.
#' @param more.cluster.args Named list of further arguments to be passed to \code{\link{clusterGraph}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#' @param output.name String containing the name of the column of the \code{\link[SummarizedExperiment]{colData}} in which to store the cluster assignments.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store extra clustering output.
#' If \code{NULL}, no extra clustering output is stored. 
#' @param graph.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store the SNN graph.
#' If \code{NULL}, the SNN graph is not stored. 
#'
#' @return \code{x} is returned with the cluster assignment for each cell stored in the \code{colData}.
#' Additional clustering output is stored in the \code{metadata}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' sce <- clusterGraph.se(sce)
#' table(sce$clusters)
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
clusterGraph.se <- function(
    x,
    num.neighbors = 10,
    num.threads = 1,
    more.build.args = list(),
    method = "multilevel",
    resolution = NULL,
    more.cluster.args = list(),
    reddim.type = "PCA",
    output.name = "clusters",
    meta.name = NULL,
    graph.name = NULL
) {
    graph.out <- .call(
        buildSnnGraph,
        list(.getTransposedReddim(x, reddim.type)),
        list(num.neighbors=num.neighbors, num.threads=num.threads, as.pointer=is.null(graph.name)),
        more.build.args
    )
    x <- .addBuildGraphResults(x, graph.out, graph.name=graph.name)

    res.args <- list()
    if (!is.null(resolution)) {
        res.args$multilevel.resolution <- resolution
        res.args$leiden.resolution <- resolution
    }

    clust.out <- .call(
        clusterGraph,
        list(graph.out),
        res.args,
        more.cluster.args
    )

    .addClusterGraphResults(x, clust.out, output.name=output.name, meta.name=meta.name)
}

.addBuildGraphResults <- function(x, graph, graph.name) {
    if (!is.null(graph.name)) {
        S4Vectors::metadata(x)[[graph.name]] <- graph
    }
    x
}

.addClusterGraphResults <- function(x, res, output.name, meta.name) {
    SummarizedExperiment::colData(x)[[output.name]] <- res$membership
    if (!is.null(meta.name)) {
        res$membership <- NULL
        S4Vectors::metadata(x)[[meta.name]] <- res
    }
    x
}
