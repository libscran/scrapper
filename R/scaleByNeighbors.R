#' Scale and combine multiple embeddings 
#'
#' Scale multiple embeddings (usually derived from different modalities for the same cells) so that their within-population variances are comparable,
#' and then combine them into a single embedding matrix for further analyses like clustering, t-SNE, etc.
#' The aim is to equalize uninteresting variance across modalities so that high technical variance in one modality does not drown out interesting biology in another modality.
#'
#' @param x List of numeric matrices of principal components or other embeddings, one for each modality.
#' For each entry, rows are dimensions and columns are cells.
#' All entries should have the same number of columns but may have different numbers of rows.
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to define the scaling factor.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} object specifying how to perform the neighbor search.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param weights Numeric vector of length equal to that of \code{x}, specifying the weights to apply to each modality.
#' Each value represents a multiplier of the within-population variance of its modality, i.e., larger values increase the contribution of that modality in the combined output matrix.
#' \code{NULL} is equivalent to an all-1 vector, i.e., all modalities are scaled to have the same within-population variance.
#'
#' @return List containing \code{scaling}, a vector of scaling factors to be aplied to each embedding;
#' and \code{combined}, a numeric matrix creating by scaling each entry of \code{x} by \code{scaling} and then \code{rbind}ing them together.
#'
#' @seealso
#' \url{https://libscran.github.io/mumosa/}, for the basis and caveats of this approach.
#'
#' @examples
#' pcs <- list(
#'     gene = matrix(rnorm(10000), ncol=200),
#'     protein = matrix(rnorm(1000, sd=3), ncol=200),
#'     guide = matrix(rnorm(2000, sd=5), ncol=200)
#' )
#'
#' out <- scaleByNeighbors(pcs)
#' out$scaling
#' dim(out$combined)
#' @author Aaron Lun
#' @export
#' @importFrom BiocNeighbors findDistance AnnoyParam
scaleByNeighbors <- function(x, num.neighbors=20, num.threads=1, weights=NULL, BNPARAM=AnnoyParam()) {
    nmod <- length(x)

    ncols <- lapply(x, ncol)
    if (length(unique(ncols)) != 1L) {
        stop("all entries of 'x' should have the same number of columns")
    }

    distances <- vector("list", nmod)
    for (i in seq_along(x)) {
        distances[[i]] <- findDistance(x[[i]], transposed=TRUE, k=num.neighbors, num.threads=num.threads, BNPARAM=BNPARAM)
    }

    scaling <- scale_by_neighbors(distances)
    if (!is.null(weights)) {
        scaling <- scaling * sqrt(weights)
    }

    names(scaling) <- names(x)
    for (i in seq_along(x)) {
        x[[i]] <- x[[i]] * scaling[i]
    }

    list(
         scaling = scaling,
         combined = do.call(rbind, x)
    )
}
