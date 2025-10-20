#' Subsample cells based on their neighbors
#'
#' Subsample a dataset by selecting cells to represent all of their nearest neighbors.
#' The aim is to preserve the relative density of the original dataset while guaranteeing representation of low-frequency subpopulations. 
#'
#' @param x A numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#'
#' Alternatively, an index constructed by \code{\link[BiocNeighbors]{buildIndex}}.
#'
#' Alternatively, a list containing existing nearest-neighbor search results.
#' This should contain:
#' \itemize{
#' \item \code{index}, an integer matrix where rows are neighbors and columns are cells.
#' Each column contains 1-based indices for the nearest neighbors of the corresponding cell, ordered by increasing distance.
#' \item \code{distance}, a numeric matrix of the same dimensions as \code{index},
#' containing the distances to each of the nearest neighbors.
#' }
#' The number of neighbors should be equal to \code{num.neighbors}, otherwise a warning is raised.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use.
#' Larger values result in stronger downsampling. 
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying the algorithm to use.
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#' @param min.remaining Integer scalar specifying the minimum number of remaining neighbors that a cell must have in order to be considered for selection.
#' This should be less than or equal to \code{num.neighbors}.
#' Larger values result in stronger downsampling. 
#' @param num.threads Integer scalar specifying the number of threads to use for the nearest-neighbor search.
#' Only used if \code{x} does not contain existing nearest-neighbor results.
#'
#' @details
#' Starting from the densest region in the high-dimensional space, we select an observation for inclusion into the subsampled dataset.
#' Every time we select an observation, we remove it and all of its nearest neighbors from the dataset.
#' We then select the next observation with the most remaining neighbors, with ties broken by density; this is repeated until there are no more observations.
#'
#' The premise is that each selected observation serves as a representative for its nearest neighbors.
#' This ensures that the subsampled points are well-distributed across the original dataset.
#' Low-frequency subpopulations will always have at least a few representatives if they are sufficiently distant from other subpopulations.
#' We also preserve the relative density of the original dataset as more representatives will be generated from high-density regions. 
#' 
#' @return Integer vector with the indices of the selected cells in the subsample.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/nenesub/}, for more details on the underlying algorithm.
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
    } else {
        .checkNeighborIndices(x$index, num.neighbors)
    }
    subsample_by_neighbors(x$index, x$distance, min_remaining=min.remaining) 
}
