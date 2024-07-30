#' Build a shared nearest neighbor graph
#'
#' Build a shared nearest neighbor (SNN) graph where each node is a cell.
#' Edges are formed between cells that share one or more nearest neighbors,
#' weighted by the number or importance of those shared neighbors.
#'
#' @param x For \code{buildSnnGraph}, a numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#' Alternatively, an index constructed by \code{\link{buildIndex}}.
#'
#' For \code{buildSnnGraphFromNeighbors}, an integer matrix where rows are neighbors and columns are cells.
#' Each column contains 1-based indices for the nearest neighbors of the corresponding cell, ordered by increasing distance.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to construct the graph.
#' @param weight.scheme String specifying the weighting scheme to use for constructing the SNN graph.
#' This can be \code{"ranked"} (default), \code{"jaccard"} or \code{"number"}.
#' @param as.pointer Logical scalar indicating whether to return a pointer to the graph.
#' This avoids an unnecessary serialization to/from R objects in downstream functions like \code{\link{clusterGraph}}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use.
#'
#' @return If \code{as.pointer=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{vertices}, an integer scalar specifying the number of vertices in the graph (i.e., cells in \code{x}).
#' \item \code{edges}, an integer vector of 1-based indices for graph edges.
#' Pairs of values represent the endpoints of an (undirected) edge,
#' i.e., \code{edges[1:2]} form the first edge, \code{edges[3:4]} form the second edge and so on.
#' \item \code{weights}, a numeric vector of weights for each edge.
#' This has length equal to half the length of \code{edges}.
#' }
#'
#' If \code{as.pointer=TRUE}, an external pointer is returned that can be used in \code{\link{clusterGraph}}. 
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{build_snn_graph} function in \url{https://libscran.github.io/scran_graph_cluster}.
#'
#' @examples
#' data <- matrix(rnorm(10000), ncol=1000)
#' out <- buildSnnGraph(data)
#' str(out)
#'
#' # We can use this to make an igraph::graph.
#' g <- igraph::make_graph(out$edges, n = out$vertices)
#' igraph::E(g)$weight <- out$weight
#'
#' @export 
#' @importFrom BiocNeighbors findKNN AnnoyParam
buildSnnGraph <- function(x, num.neighbors=10, weight.scheme="ranked", as.pointer=FALSE, num.threads=1, BNPARAM=AnnoyParam()) {
    neighbors <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance=FALSE, num.threads=num.threads, BNPARAM=BNPARAM)
    buildSnnGraphFromNeighbors(neighbors$index, weight.scheme=weight.scheme, as.pointer=as.pointer, num.threads=num.threads)
}

#' @export
#' @rdname buildSnnGraph
buildSnnGraphFromNeighbors <- function(x, weight.scheme="ranked", as.pointer=FALSE, num.threads=1) {
    build_snn_graph(x, scheme=weight.scheme, num_threads=num.threads, raw=as.pointer)
}
