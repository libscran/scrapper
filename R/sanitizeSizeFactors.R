#' Sanitize size factors
#'
#' Replace invalid size factors, i.e., zero, negative, infinite or NaN values.
#'
#' @param size.factors Numeric vector of size factors across cells.
#' @param replace.zero Logical scalar indicating whether to replace size factors of zero with the lowest positive factor.
#' If \code{FALSE}, zeros are retained.
#' @param replace.negative Logical scalar indicating whether to replace negative size factors with the lowest positive factor.
#' If \code{FALSE}, negative values are retained.
#' @param replace.infinite Logical scalar indicating whether to replace infinite size factors with the largest positive factor.
#' If \code{FALSE}, infinite values are retained.
#' @param replace.nan Logical scalar indicating whether to replace NaN size factors with unity.
#' If \code{FALSE}, NaN values are retained.
#' 
#' @return Numeric vector of length equal to \code{size.factors}, containing the sanitized size factors.
#'
#' @author Aaron Lun
#'
#' @examples
#' sf <- 2^rnorm(100)
#' sf[1] <- 0
#' sf[2] <- -1
#' sf[3] <- Inf
#' sf[4] <- NaN
#' sanitizeSizeFactors(sf)
#'
#' @seealso
#' The \code{sanitize_size_factors} function in \url{https://libscran.github.io/scran_norm/}, for more details on the sanitization.
#'
#' @export
sanitizeSizeFactors <- function(size.factors, replace.zero=TRUE, replace.negative=TRUE, replace.infinite=TRUE, replace.nan=TRUE) {
    sanitize_size_factors(size.factors, handle_zero=replace.zero, handle_negative=replace.negative, handle_infinite=replace.infinite, handle_nan=replace.nan)
}
