#' Sanitize size factors
#'
#' Replace invalid size factors, i.e., zero, negative, infinite or NaN values.
#' Such size factors can occasionally arise if, e.g., insufficient quality control was performed upstream.
#' Removing them ensures that the normalized values from \code{\link{normalizeCounts}} remain finite for sensible downstream processing.
#'
#' @param size.factors Numeric vector of size factors across cells.
#' @param replace.zero Logical scalar indicating whether to replace size factors of zero with the lowest positive factor in \code{size.factors}.
#' This ensures that the normalized values will be large to reflect the extremity of the scaling, but still finite for sensible downstream processing. 
#' If \code{FALSE}, zeros are retained.
#' @param replace.negative Logical scalar indicating whether to replace negative size factors with the lowest positive factor in \code{size.factors}.
#' This ensures that the normalized values will be large to reflect the extremity of the scaling, but still finite for sensible downstream processing. 
#' If \code{FALSE}, negative values are retained.
#' @param replace.infinite Logical scalar indicating whether to replace infinite size factors with the largest positive factor in \code{size.factors}.
#' This ensures that any normalized values will be, at least, finite; the choice of a relatively large replacement value reflects the extremity of the scaling.
#' If \code{FALSE}, infinite values are retained.
#' @param replace.nan Logical scalar indicating whether to replace NaN size factors with unity, e.g., scaling normalization is a no-op.
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
#' The \code{sanitize_size_factors} function in \url{https://libscran.github.io/scran_norm/}. 
#'
#' @export
sanitizeSizeFactors <- function(size.factors, replace.zero=TRUE, replace.negative=TRUE, replace.infinite=TRUE, replace.nan=TRUE) {
    sanitize_size_factors(size.factors, handle_zero=replace.zero, handle_negative=replace.negative, handle_infinite=replace.infinite, handle_nan=replace.nan)
}
