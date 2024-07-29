#' Build a shared nearest neighbor graph
#'
#' Build a shared nearest neighbor (SNN) graph where each node is a cell.
#' Edges are formed between cells that share one or more nearest neighbors,
#' weighted by the number or importance of those shared neighbors.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#' Alternatively, an index constructed by \code{\link{buildIndex}}.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to construct the graph.
#' @param weight.scheme String specifying the weighting scheme to use for constructing the SNN graph.
#' This can be \code{"ranked"} (default), \code{"jaccard"} or \code{"number"}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use.
#'
#' @return List containing:
#' \itemize{
#' \item \code{edges}, an integer vector of 0-based indices for graph edges.
#' Pairs of values represent the endpoints of an (undirected) edge,
#' i.e., \code{edges[1:2]} form the first edge, \code{edges[3:4]} form the second edge and so on.
#' \item \code{weights}, a numeric vector of weights for each edge.
#' This has length equal to half the length of \code{edges}.
#' }
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
#' g <- igraph::make_graph(out$edges + 1L, weights=out$weights)
#'
#' @export 
#' @importFrom BiocNeighbors findKNN AnnoyParam
buildSnnGraph <- function(x, num.neighbors=10, weight.scheme="ranked", num.threads=1, BNPARAM=AnnoyParam()) {
    neighbors <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance=FALSE, num.threads=num.threads, BNPARAM=BNPARAM)
    build_snn_graph(neighbors$index - 1L, scheme=weight.scheme, num_threads=num.threads)
}
