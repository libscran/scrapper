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
#'
#' If \code{NULL}, the default value in \code{\link{normalizeCountsDefaults}} is used.
#' @param pseudo.count Numeric scalar specifying the positive pseudo-count to add before any log-transformation.
#' Larger values shrink the differences between cells towards zero, reducing spurious differences (but also signal) at low counts - see \code{\link{choosePseudoCount}} for comments.
#'
#' If \code{NULL}, the default value in \code{\link{normalizeCountsDefaults}} is used.
#'
#' This argument is ignored if \code{log = FALSE}.
#' @param log.base Numeric scalar specifying the base of the log-transformation.
#' This is typically set to 2 or 10 for convenient interpretation of the log-values.
#'
#' If \code{NULL}, the default value in \code{\link{normalizeCountsDefaults}} is used.
#'
#' This argument is ignored if \code{log = FALSE}.
#' @param preserve.sparsity Logical scalar indicating whether to preserve sparsity for \code{pseudo.count != 1}.
#' This improves memory and compute efficiency in downstream steps but changes the interpretation of the returned matrix. 
#' Specifically, if \code{TRUE}, users should manually add \code{log(pseudo.count, log.base)} to the returned matrix to obtain the desired log-transformed expression values.
#'
#' If \code{NULL}, the default value in \code{\link{normalizeCountsDefaults}} is used.
#'
#' This argument is ignored if \code{log = FALSE} or \code{pseudo.count = 1}.
#'
#' @return A normalized expression matrix in varying forms:
#' \itemize{
#' \item If \code{x} is a matrix-like object, a matrix-like object is returned containing the (log-transformed) normalized expression matrix.
#' If \code{delayed=TRUE}, the returned object is a \link[DelayedArray]{DelayedArray}, otherwise the return type depends on the type of \code{x} and the operations involved.
#' \item If \code{x} is an external pointer produced by \code{\link[beachmat]{initializeCpp}}, a new external pointer is returned containing the normalized expression matrix.
#' }
#'
#' @seealso
#' The \code{normalize_counts} function in \url{https://libscran.github.io/scran_norm/}.
#'
#' \code{\link{normalizeRnaCounts.se}} and related functions, which compute normalized expression values from a \link[SummarizedExperiment]{SummarizedExperiment}.
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
normalizeCounts <- function(x, size.factors, log = NULL, pseudo.count = NULL, log.base = NULL, preserve.sparsity = NULL, delayed=TRUE) {
    .checkSEX(x, "normalizeRnaCounts.se")
    if (is(x, "externalptr")) {
        return(normalize_counts(x, size.factors, log=log, log_base=log.base, pseudo_count=pseudo.count, preserve_sparsity=preserve.sparsity))
    }

    def <- normalizeCountsDefaults()
    if (is.null(log)) {
        log <- def$log
    }
    if (is.null(pseudo.count)) {
        pseudo.count <- def$pseudo.count
    }
    if (is.null(log.base)) {
        log.base <- def$log.base
    }
    if (is.null(preserve.sparsity)) {
        preserve.sparsity <- def$preserve.sparsity
    }

    if (!log) {
        if (delayed) {
            x <- DelayedArray(x)
        }
        return(t(t(x)/size.factors))
    }

    if (pseudo.count != 1 && preserve.sparsity) {
        size.factors <- size.factors * pseudo.count
        pseudo.count <- 1
    }

    if (!delayed) {
        .log_normalize(x, size.factors, pseudo.count=pseudo.count, log.base=log.base)
    } else {
        LogNormalizedMatrix(x, size.factors, pseudo.count=pseudo.count, log.base=log.base)
    }
}

#' Default parameters for \code{\link{normalizeCounts}}
#' @description Default parameters from the underlying C++ library.
#' These may be overridden by defaults in the \code{\link{normalizeCounts}} function signature.
#' @return Named list of default values for various function arguments. 
#' @author Aaron Lun
#' @examples
#' normalizeCountsDefaults()
#' @export
normalizeCountsDefaults <- function() normalize_counts_defaults()
