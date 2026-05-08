#' Compute size factors for ADT counts
#'
#' Compute size factors from an ADT count matrix using the CLRm1 method.
#' This is a variant of the centered log-ratio (CLR) method, where the size factors are defined from the geometric mean of counts within each cell.
#'
#' @param x A matrix-like object containing ADT count data.
#' Rows correspond to tags and columns correspond to cells.
#' @param num.threads Number of threads to use.
#'
#' If \code{NULL}, the default value in \code{\link{computeClrm1FactorsDefaults}} is used.
#'
#' @return Numeric vector containing the CLRm1 size factor for each cell.
#' Note that these size factors are not centered and should be passed through, e.g., \code{\link{centerSizeFactors}} before normalization.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://github.com/libscran/clrm1}, for a description of the CLRm1 method.
#'
#' \code{\link{normalizeAdtCounts.se}}, which computes CLRm1 factors prior to normalization.
#'
#' @examples
#' library(Matrix)
#' x <- abs(rsparsematrix(1000, 100, 0.1) * 10)
#' head(computeClrm1Factors(x))
#'
#' @export
#' @importFrom beachmat initializeCpp
computeClrm1Factors <- function(x, num.threads = NULL) {
    .checkSEX(x, "normalizeAdtCounts.se")
    ptr <- initializeCpp(x, .check.na=FALSE)
    compute_clrm1_factors(ptr, num_threads=num.threads)
}

#' Default parameters for \code{\link{computeClrm1Factors}}
#' @description Default parameters from the underlying C++ library.
#' These may be overridden by defaults in the \code{\link{computeClrm1Factors}} function signature.
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' computeClrm1FactorsDefaults()
#' @export
computeClrm1FactorsDefaults <- function() compute_clrm1_factors_defaults()
