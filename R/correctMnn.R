#' Batch correction with mutual nearest neighbors
#'
#' Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing low-dimensional coordinates (e.g., from \code{\link{runPca}}).
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use when identifying MNN pairs.
#' @param num.mads Numeric scalar specifying the number of median absolute deviations to use for removing outliers in the center-of-mass calculations.
#' @param robust.iterations Integer scalar specifying the number of iterations for robust calculation of the center of mass.
#' @param robust.trim Numeric scalar in [0, 1) specifying the trimming proportion for robust calculation of the center of mass.
#' @param mass.cap Integer scalar specifying the cap on the number of observations to use for center-of-mass calculations on the reference dataset.
#' A value of 100,000 may be appropriate for speeding up correction of very large datasets.
#' If \code{NULL}, no cap is used.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param order Vector containing levels of \code{batch} in the desired merge order.
#' If \code{NULL}, a suitable merge order is automatically determined.
#' @param reference.policy String specifying the policy to use to choose the first reference batch.
#' This can be based on the largest batch (\code{"max-size"}),
#' the most variable batch (\code{"max-variance"}), 
#' the batch with the largest residual sum of squares (\code{"max-rss"}),
#' or the first specified input (\code{"input"}).
#' Only used for automatic merges, i.e., when \code{order=NULL}. 
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying the nearest-neighbor algorithm to use.
#'
#' @return List containing:
#' \itemize{
#' \item \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
#' \item \code{merge.order}, character vector containing the unique levels of \code{batch} in the automatically determined merge order.
#' The first level in this vector is used as the reference batch; all other batches are iteratively merged to it.
#' \item \code{num.pairs}, integer vector of length equal to the number of batches minus 1.
#' This contains the number of MNN pairs at each merge.
#' }
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
    num.mads=3,
    robust.iterations=2,
    robust.trim=0.25,
    mass.cap=NULL,
    order=NULL,
    reference.policy=c("max-rss", "max-size", "max-variance", "input"),
    BNPARAM=AnnoyParam(),
    num.threads=1)
{
    block <- .transformFactor(block)

    if (!is.null(order)) {
        order <- match(order, block$names)
        if (!identical(sort(order), seq_along(order))) {
            stop("'order' should contain unique values in 'block'"); 
        }
        order <- order - 1L
    }

    if (is.null(mass.cap)) {
        mass.cap <- -1
    }

    output <- correct_mnn(
        x, 
        block$index, 
        num_neighbors=num.neighbors,
        num_mads=num.mads,
        robust_iterations=robust.iterations,
        robust_trim=robust.trim,
        mass_cap=mass.cap,
        num_threads=num.threads, 
        order=order, 
        ref_policy=match.arg(reference.policy), 
        builder=defineBuilder(BNPARAM)$builder
    )

    output$merge.order <- block$names[output$merge.order + 1L]
    output
}
