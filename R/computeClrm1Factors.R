#' Compute size factors for ADT counts
#'
#' Compute size factors from an ADT count matrix using the CLRm1 method.
#'
#' @param x A matrix-like object containing ADT count data.
#' Rows correspond to tags and columns correspond to cells.
#' @param num.threads Number of threads to use.
#'
#' @return Numeric vector containing the CLRm1 size factor for each cell.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://github.com/libscran/clrm1}, for a description of the CLRm1 method.
#'
#' @examples
#' library(Matrix)
#' x <- abs(rsparsematrix(1000, 100, 0.1) * 10)
#' head(computeClrm1Factors(x))
#'
#' @export
#' @importFrom beachmat initializeCpp
computeClrm1Factors <- function(x, num.threads=1) {
    ptr <- initializeCpp(x, .check.na=FALSE)
    compute_clrm1_factors(ptr, num_threads=num.threads)
}
