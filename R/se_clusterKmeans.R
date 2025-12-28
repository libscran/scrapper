#' k-means clustering of cells in a SingleCellExperiment
#'
#' Perform k-means clustering on an existing low-dimensional embedding
#' by calling \code{\link{clusterKmeans}} on a reduced dimension entry in a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param k Number of clusters, passed to \code{\link{clusterKmeans}}.
#' @param num.threads Number of threads, passed to \code{\link{clusterKmeans}}.
#' @param more.kmeans.args Named list of further arguments to be passed to \code{\link{clusterKmeans}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' Alternatively, a named integer or character vector of length 1,
#' where the name specifies an alternative experiment of \code{x} and the value is the name/index of a \code{reducedDim} entry in that alternative experiment.
#' @param output.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the cluster assignments.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store extra clustering output.
#' If \code{NULL}, no extra clustering output is stored. 
#'
#' @return \code{x} is returned with the cluster assignment for each cell stored in the \code{colData}.
#' Additional clustering output is stored in the \code{metadata}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' sce <- clusterKmeans.se(sce, k=10)
#' table(sce$clusters)
#' 
#' @export
#' @importFrom S4Vectors metadata metadata<-
clusterKmeans.se <- function(x, k, num.threads=1, more.kmeans.args = list(), reddim.type = "PCA", output.name = "clusters", meta.name = NULL) {
    clout <- .call(
        clusterKmeans,
        list(.getTransposedReddim(x, reddim.type)),
        list(k=k, num.threads=num.threads),
        more.kmeans.args
    )

    SummarizedExperiment::colData(x)[[output.name]] <- clout$clusters
    if (!is.null(meta.name)) {
        clout$clusters <- NULL
        metadata(x)[[meta.name]] <- clout
    }

    x
}
