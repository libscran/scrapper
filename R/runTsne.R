#' t-stochastic neighbor embedding
#'
#' Compute t-SNE coordinates to visualize similarities between cells.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#'
#' Alternatively, a named list of nearest-neighbor search results.
#' This should contain \code{index}, an integer matrix where rows are neighbors and columns are cells.
#' Each column contains 1-based indices for the nearest neighbors of the corresponding cell, ordered by increasing distance.
#' The number of neighbors should be the same as \code{tsnePerplexityToNeighbors(perplexity)}.
#'
#' Alternatively, an index constructed by \code{\link{buildIndex}}.
#' @param perplexity Numeric scalar specifying the perplexity to use in the t-SNE algorithm.
#' @param max.depth Integer scalar specifying the maximum depth of the Barnes-Hut quadtree.
#' Smaller values (7-10) improve speed at the cost of accuracy.
#' @param leaf.approximation Logical scalar indicating whether to use the \dQuote{leaf approximation} approach,
#' which sacrifices some accuracy for greater speed.
#' Only effective when \code{max.depth} is small enough for multiple cells to be assigned to the same leaf node of the quadtree.
#' @param seed Integer scalar specifying the seed to use for generating the initial coordinates.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param max.iterations Integer scalar specifying the maximum number of iterations to perform.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the algorithm to use.
#' Only used if \code{x} is not a list of existing nearest-neighbor search results.
#' 
#' @return 
#' For \code{runTsne}, a numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#'
#' For \code{tsnePerplexityToNeighbors}, an integer scalar specifying the number of neighbors to use for a given \code{perplexity}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/qdtsne/}, for an explanation of the approximations.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runTsne(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
runTsne <- function(x, perplexity=30, max.depth=20, leaf.approximation=FALSE, max.iterations=500, seed=42, num.threads=1, BNPARAM=AnnoyParam()) {
    if (!is.list(x)) {
        k <- perplexity_to_neighbors(perplexity)
        x <- findKNN(x, k=k, transposed=TRUE, get.index="transposed", get.distance="transposed", num.threads=num.threads, BNPARAM=BNPARAM)
    } else {
        .checkIndices(x$index)
    }

    output <- run_tsne(
        nnidx=x$index, 
        nndist=x$distance, 
        perplexity=perplexity,
        leaf_approx=leaf.approximation,
        max_depth=max.depth,
        max_iter=max.iterations,
        num_threads=num.threads,
        seed=seed
    )

    t(output)
}

#' @export
#' @rdname runTsne
tsnePerplexityToNeighbors <- function(perplexity) {
    perplexity_to_neighbors(perplexity)
}
