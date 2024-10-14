#' Uniform manifold approxation and projection
#'
#' Compute UMAP coordinates to visualize similarities between cells.
#'
#' @inheritParams runTsne
#' @param num.dim Integer scalar specifying the number of dimensions of the output embedding.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use in the UMAP algorithm.
#' @param min.dist Numeric scalar specifying the minimum distance between points.
#' @param seed Integer scalar specifying the seed to use. 
#' @param num.epochs Integer scalar specifying the number of epochs to perform.
#' If set to -1, an appropriate number of epochs is chosen based on \code{ncol(x)}.
#' @param parallel.optimization Logical scalar specifying whether to parallelize the optimization step.
#' 
#' @return 
#' A numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/umappp/}, for details on the underlying implementation.
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runUmap(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
runUmap <- function(x, 
    num.dim=2,
    num.neighbors=15, 
    num.epochs=-1, 
    min.dist=0.1, 
    seed=1234567890, 
    num.threads=1,
    parallel.optimization=FALSE,
    BNPARAM=AnnoyParam())
{
    if (!is.list(x)) {
        x <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance="transposed", num.threads=num.threads, BNPARAM=BNPARAM)
    } else {
        .checkIndices(x$index, num.neighbors)
    }

    output <- run_umap(
        nnidx=x$index, 
        nndist=x$distance, 
        ndim=num.dim,
        min_dist=min.dist,
        num_epochs=num.epochs,
        num_threads=num.threads,
        parallel_optimization=parallel.optimization,
        seed=seed
    )

    t(output)
}
