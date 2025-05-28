#' Batch correction with mutual nearest neighbors
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing low-dimensional coordinates (e.g., from \code{\link{runPca}}).
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use when identifying MNN pairs.
#' @param num.steps Integer scalar specifying the number of steps for the recursive neighbor search to compute the center of mass.
#' @param num.mads Deprecated and ignored.
#' @param robust.iterations Deprecated and ignored.
#' @param robust.trim Deprecated and ignored.
#' @param mass.cap Deprecated and ignored.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param order Deprecated and ignored, the merge order is now always automatically determined.
#' @param reference.policy Deprecated, use \code{merge.policy} instead. 
#' @param merge.policy String specifying the policy to use to choose the order of batches to merge.
#' This can be based on the size of the batch (\code{"size"}),
#' the variance within each batch (\code{"variance"}), 
#' the residual sum of squares of each batch (\code{"rss"}),
#' or the input order (\code{"input"}).
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying the nearest-neighbor algorithm to use.
#'
#' @return List containing \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
#'
#' @seealso
#' \url{https://libscran.github.io/mnncorrect/}, for more details on the underlying implementation.
#'
#' @author Aaron Lun
#'
#' @examples
#' # Mocking up a dataset with multiple batches.
#' x <- matrix(rnorm(10000), nrow=10)
#' b <- sample(3, ncol(x), replace=TRUE)
#' x[,b==2] <- x[,b==2] + 3
#' x[,b==3] <- x[,b==3] + 5
#' lapply(split(colMeans(x), b), mean) # indeed the means differ...
#'
#' corrected <- correctMnn(x, b)
#' str(corrected)
#' lapply(split(colMeans(corrected$corrected), b), mean) # now merged.
#'
#' @export
#' @importFrom BiocNeighbors defineBuilder
correctMnn <- function(
    x,
    block,
    num.neighbors=15,
    num.steps=1,
    merge.policy=c("rss", "size", "variance", "input"),
    num.mads=NULL,
    robust.iterations=NULL,
    robust.trim=NULL,
    mass.cap=NULL,
    order=NULL,
    reference.policy=NULL,
    BNPARAM=AnnoyParam(),
    num.threads=1)
{
    block <- .transformFactor(block)

    if (!is.null(num.mads)) {
        .Deprecated(old="num.mads=")
    }
    if (!is.null(robust.iterations)) {
        .Deprecated(old="robust.iterations=")
    }
    if (!is.null(robust.trim)) {
        .Deprecated(old="robust.trim=")
    }
    if (!is.null(mass.cap)) {
        .Deprecated(old="mass.cap=")
    }
    if (!is.null(order)) {
        .Deprecated(old="order=")
    }

    if (!is.null(reference.policy)) {
        merge.policy <- sub("^max-", "", reference.policy)
        .Deprecated(old=sprintf("reference.policy=%s", deparse(reference.policy)), new=sprintf("merge.policy=%s", deparse(merge.policy)))
    }
    merge.policy <- match.arg(merge.policy)

    output <- correct_mnn(
        x, 
        block$index, 
        num_neighbors=num.neighbors,
        num_steps=num.steps,
        num_threads=num.threads,
        merge_policy=merge.policy,
        builder=defineBuilder(BNPARAM)$builder
    )

    output
}
