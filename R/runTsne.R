#' t-stochastic neighbor embedding
#'
#' Compute t-SNE coordinates to visualize similarities between cells.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#'
#' Alternatively, a named list of nearest-neighbor search results like that returned by \code{\link[BiocNeighbors]{findKNN}}.
#' This should contain \code{index}, an integer matrix where rows are neighbors and columns are cells;
#' and \code{distance}, a numeric matrix of the same dimensions containing the distances to each neighbor.
#' Each column contains 1-based indices for the nearest neighbors of the corresponding cell, ordered by increasing distance.
#' The number of neighbors should be the same as \code{num.neighbors}, otherwise a warning is raised.
#'
#' Alternatively, an index constructed by \code{\link[BiocNeighbors]{buildIndex}}.
#' @param perplexity Numeric scalar specifying the perplexity to use in the t-SNE algorithm.
#' Higher perplexities will focus on global structure, at the cost of increased runtime and decreased local resolution.
#' @param num.neighbors Integer scalar specifying the number of neighbors, typically derived from \code{perplexity}.
#' If \code{x} contains pre-computed neighbor search results with a different number of neighbors than \code{num.neighbors}, an error is thrown;
#' this can be suppressed by setting \code{num.neighbors = NULL}.
#' @param theta Numeric scalar specifying the approximation level for the Barnes-Hut calculation of repulsive forces.
#' Lower values increase accuracy at the cost of increased compute time.
#' All values should be non-negative.
#' @param early.exaggeration.iterations Integer scalar specifying the number of iterations of the early exaggeration phase,
#' where clusters are artificially compacted to leave more empty space so that cells can easily relocate to find a good global organization.
#' Larger values improve convergence within this phase at the cost of reducing the remaining iterations in \code{max.iterations}.
#' @param exaggeration.factor Numeric scalar containing the exaggeration factor for the early exaggeration phase (see \code{early.exaggeration.iterations}).
#' Larger values increase the attraction between nearest neighbors to favor local structure.
#' @param momentum.switch.iterations Integer scalar specifying the number of iterations to perform before switching from the starting momentum to the final momentum.
#' Higher momentums can improve convergence by increasing the step size and smoothing over local oscillations, at the risk of potentially skipping over relevant minima.
#' @param start.momentum Numeric scalar containing the starting momentum, to be used in the iterations before the momentum switch at \code{momentum.switch.iterations}.
#' This is usually lower than \code{final.momentum} to avoid skipping over suitable local minima. 
#' @param final.momentum Numeric scalar containing the final momentum, to be used in the iterations after the momentum switch at \code{momentum.switch.iterations}.
#' This is usually higher than \code{start.momentum} to accelerate convergence to the local minima once the observations are moderately well-organized.
#' @param eta Numeric scalar containing the learning rate, used to scale the updates for each cell.
#' Larger values can speed up convergence at the cost of skipping over local minima. 
#' @param max.depth Integer scalar specifying the maximum depth of the Barnes-Hut quadtree.
#' If neighboring cells cannot be separated before the maximum depth is reached, they will be assigned to the same leaf node of the quadtree.
#' Smaller values (7-10) improve speed by bounding the recursion depth at the cost of accuracy.
#' @param leaf.approximation Logical scalar indicating whether to use the \dQuote{leaf approximation}.
#' If \code{TRUE}, repulsive forces are computed between leaf nodes and re-used across all cells assigned to that leaf node. 
#' This sacrifices some accuracy for greater speed, assuming that \code{max.depth} is small enough for multiple cells to be assigned to the same leaf.
#' @param seed Integer scalar specifying the seed to use for generating the initial coordinates.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param max.iterations Integer scalar specifying the maximum number of iterations to perform.
#' Larger values improve convergence at the cost of compute time.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying the algorithm to use.
#' Only used if \code{x} is not a prebuilt index or a list of existing nearest-neighbor search results.
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
#' \code{\link{runTsne.se}}, to run t-SNE on a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runTsne(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @references
#' van der Maaten LJP and Hinton GE (2008). 
#' Visualizing high-dimensional data using t-SNE. 
#' \emph{Journal of Machine Learning Research_} 9, 2579-2605.
#' 
#' van der Maaten LJP (2014). 
#' Accelerating t-SNE using tree-based algorithms. 
#' \emph{Journal of Machine Learning Research} 15, 3221-3245.
#' 
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
runTsne <- function(
    x,
    perplexity=30,
    num.neighbors=tsnePerplexityToNeighbors(perplexity),
    theta=1,
    early.exaggeration.iterations=250,
    exaggeration.factor=12,
    momentum.switch.iterations=250,
    start.momentum=0.5,
    final.momentum=0.8,
    eta=200,
    max.depth=7,
    leaf.approximation=FALSE,
    max.iterations=500,
    seed=42,
    num.threads=1,
    BNPARAM=AnnoyParam())
{
    .checkSEX(x, "runTsne.se")

    if (!is.list(x)) {
        x <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance="transposed", num.threads=num.threads, BNPARAM=BNPARAM)
    } else {
        .checkNeighborResults(x$index, x$distance)
        if (!is.null(num.neighbors) && nrow(x$index) != num.neighbors) {
            warning("'nrow(x$index)' is not consistent with 'num.neighbors'")
        }
    }

    output <- run_tsne(
        nnidx=x$index, 
        nndist=x$distance, 
        perplexity=perplexity,
        theta=theta,
        early_exaggeration_iterations=early.exaggeration.iterations,
        exaggeration_factor=exaggeration.factor,
        momentum_switch_iterations=momentum.switch.iterations,
        start_momentum=start.momentum,
        final_momentum=final.momentum,
        eta=eta,
        max_depth=max.depth,
        leaf_approx=leaf.approximation,
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
