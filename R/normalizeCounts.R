#' Normalize the count matrix
#'
#' Apply scaling normalization and log-transformation to a count matrix.
#' This yields normalized expression values that can be used in downstream procedures like PCA. 
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Values are expected to be non-negative counts.
#' Alternatively, an external pointer created by \code{\link[beachmat]{initializeCpp}}.
#' @param size.factors A numeric vector of length equal to the number of cells in \code{x}, containing positive size factors for all cells.
#' Any invalid values should be replaced with \code{\link{sanitizeSizeFactors}}.
#' For most applications, these size factors should also be centered with \code{\link{centerSizeFactors}}.
#' @param delayed Logical scalar indicating whether operations on a matrix-like \code{x} should be delayed.
#' This improves memory efficiency at the cost of some speed in downstream operations.
#' @param log Logical scalar indicating whether log-transformation should be performed.
#' This ensures that downstream analyses like t-tests and distance calculations focus on relative fold-changes rather than absolute differences.
#' The log-transformation also provides some measure of variance stabilization so that the downstream analyses are not dominated by sampling noise at large counts.
#' @param pseudo.count Numeric scalar specifying the positive pseudo-count to add before any log-transformation.
#' Larger values shrink the differences between cells towards zero, reducing spurious differences (but also signal) at low counts - see \code{\link{choosePseudoCount}} for comments.
#' Ignored if \code{log=FALSE}.
#' @param log.base Numeric scalar specifying the base of the log-transformation.
#' Ignored if \code{log=FALSE}.
#' @param preserve.sparsity Logical scalar indicating whether to preserve sparsity for \code{pseudo.count != 1}.
#' If \code{TRUE}, users should manually add \code{log(pseudo.count, log.base)} to the returned matrix to obtain the desired log-transformed expression values.
#' Ignored if \code{log = FALSE} or \code{pseudo.count = 1}.
#'
#' @return If \code{x} is a matrix-like object and \code{delayed=TRUE}, a \link[DelayedArray]{DelayedArray} is returned containing the (log-transformed) normalized expression matrix.
#' If \code{delayed=FALSE}, the type of the (log-)normalized matrix will depend on the operations applied to \code{x}.
#'
#' If \code{x} is an external pointer produced by \code{\link[beachmat]{initializeCpp}}, a new external pointer is returned containing the normalized expression matrix.
#'
#' @seealso
#' The \code{normalize_counts} function in \url{https://libscran.github.io/scran_norm/}.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' sf <- centerSizeFactors(colSums(x))
#' normed <- normalizeCounts(x, size.factors=sf)
#' normed
#'
#' # Passing a pointer.
#' ptr <- beachmat::initializeCpp(x)
#' optr <- normalizeCounts(ptr, sf)
#' optr
#'
#' @export
#' @importFrom methods is
#' @importFrom DelayedArray DelayedArray t
normalizeCounts <- function(x, size.factors, log=TRUE, pseudo.count=1, log.base=2, preserve.sparsity=FALSE, delayed=TRUE) {
    if (is(x, "externalptr")) {
        return(normalize_counts(x, size.factors, log=log, log_base=log.base, pseudo_count=pseudo.count, preserve_sparsity=preserve.sparsity))
    }

    if (log && pseudo.count != 1 && preserve.sparsity) {
        size.factors <- size.factors * pseudo.count
        pseudo.count <- 1
    }

    if (delayed) {
        x <- DelayedArray(x)
    }

    normalized <- t(t(x) / size.factors)
    if (!log) {
        return(normalized)
    }

    if (pseudo.count == 1) {
        return(log1p(normalized) / log(log.base))
    }

    log(normalized + pseudo.count, log.base)
}
