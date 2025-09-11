#' Batch correction with mutual nearest neighbors
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional embedding.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing coordinates in a low-dimensional embedding (e.g., from \code{\link{runPca}}).
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' @param num.neighbors Integer scalar specifying the number of neighbors in the various search steps.
#' Larger values improve the stability of the correction by increasing the number of MNN pairs and including more observations in each center of mass.
#' However, this comes at the cost of reduced resolution when matching subpopulations across batches.
#' @param num.steps Integer scalar specifying the number of steps for the recursive neighbor search to compute the center of mass.
#' Larger values mitigate the kissing effect but increase the risk of including inappropriately distant subpopulations into the center of mass.
#' @param num.mads Deprecated and ignored.
#' @param robust.iterations Deprecated and ignored.
#' @param robust.trim Deprecated and ignored.
#' @param mass.cap Deprecated and ignored.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param order Deprecated and ignored, the merge order is now always automatically determined.
#' @param reference.policy Deprecated, use \code{merge.policy} instead. 
#' @param merge.policy String specifying the policy to use to choose the order of batches to merge.
#' \itemize{
#' \item \code{"input"} will use the input order of the batches.
#' Observations in the last batch are corrected first, and then the second-last batch, and so on.
#' This allows users to control the merge order by simply changing the inputs.
#' \item \code{"size"} will merge batches in order of increasing size (i.e., the number of observations).
#' So, the smallest batch is corrected first while the largest batch is unchanged.
#' The aim is to lower compute time by reducing the number of observations that need to be reprocessed in later merge steps.
#' \item \code{"variance"} will merge batches in order of increasing variance between observations. 
#' So, the batch with the lowest variance is corrected first while the batch with the highest variance is unchanged.
#' The aim is to lower compute time by encouraging more observations to be corrected to the most variable batch, thus avoid reprocessing in later merge steps.
#' \item \code{"rss"} will merge batches in order of increasing residual sum of squares (RSS).
#' This is effectively a compromise between \code{"variance"} and \code{"size"}.
#' }
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying the nearest-neighbor algorithm to use.
#'
#' @return List containing \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
#'
#' @seealso
#' \url{https://libscran.github.io/mnncorrect/}, for more details on the underlying implementation.
#'
#' @references 
#' Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018).
#' Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.
#' \emph{Nat. Biotechnol.} 36(5):421-427
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
