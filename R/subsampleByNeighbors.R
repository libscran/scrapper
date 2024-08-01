#' Subsample cells based on their neighbors
#'
#' Subsample a dataset by selecting cells to represent all of their nearest neighbors.
#'
#' @param x A numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#'
#' Alternatively, an index constructed by \code{\link{buildIndex}}.
#'
#' Alternatively, a list containing existing nearest-neighbor search results.
#' This should contain:
#' \itemize{
#' \item \code{index}, an integer matrix where rows are neighbors and columns are cells.
#' Each column contains 1-based indices for the nearest neighbors of the corresponding cell, ordered by increasing distance.
#' \item \code{distance}, a numeric matrix of the same dimensions as \code{index},
#' containing the distances to each of the nearest neighbors.
#' }
#' @param num.neighbors Integer scalar specifying the number of neighbors to use.
#' Larger values result in greater downsampling. 
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use.
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#' @param min.remaining Integer scalar specifying the minimum number of remaining (i.e., unselected) neighbors that a cell must have in order to be considered for selection.
#' This should be less than or equal to \code{num.neighbors}.
#' @param num.threads Integer scalar specifying the number of threads to use for the nearest-neighbor search.
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#' 
#' @return Integer vector with the indices of the selected cells in the subsample.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{compute} function in \url{https://libscran.github.io/nenesub}.
#'
#' @examples
#' x <- matrix(rnorm(10000), nrow=2)
#' keep <- subsampleByNeighbors(x, 10)
#' plot(x[1,], x[2,])
#' points(x[1,keep], x[2,keep], col="red")
#' legend('topright', col=c('black', 'red'), legend=c('all', 'subsample'), pch=1)
#'
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
subsampleByNeighbors <- function(x, num.neighbors=20, min.remaining=10, num.threads=1, BNPARAM=AnnoyParam()) {
    if (!is.list(x)) {
        x <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance="transposed", num.threads=num.threads, BNPARAM=BNPARAM)
    }
    subsample_by_neighbors(x$index, x$distance, min_remaining=min.remaining) 
}
