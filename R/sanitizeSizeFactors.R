#' Sanitize size factors
#'
#' Replace invalid size factors, i.e., zero, negative, infinite or NaN values.
#' Such size factors can occasionally arise if, e.g., insufficient quality control was performed upstream.
#' Removing them ensures that the normalized values from \code{\link{normalizeCounts}} remain finite for sensible downstream processing.
#'
#' @param size.factors Numeric vector of size factors across cells.
#' @param replace.zero Deprecated, use \code{handle.zero} instead.
#' @param handle.zero String specifying how to handle replace size factors of zero.
#' If \code{"sanitize"}, zero size factors are replaced with the lowest positive factor in \code{size.factors}.
#' This ensures that the normalized values will be large to reflect the extremity of the scaling, but still finite for sensible downstream processing. 
#' If \code{"error"}, an error is raised upon encountering a size factor of zero.
#' If \code{"ignore"}, no action is taken.
#' 
#' If \code{NULL}, the default value in \code{\link{sanitizeSizeFactorsDefaults}} is used.
#' @param replace.negative Deprecated, use \code{handle.negative} instead.
#' @param handle.negative String specifying how to handle negative size factors.
#' If \code{"sanitize"}, negative size factors are replaced with the lowest positive factor in \code{size.factors}.
#' This ensures that the normalized values will be large to reflect the extremity of the scaling, but still finite for sensible downstream processing. 
#' If \code{"error"}, an error is raised upon encountering a negative size factor. 
#' If \code{"ignore"}, no action is taken.
#' 
#' If \code{NULL}, the default value in \code{\link{sanitizeSizeFactorsDefaults}} is used.
#' @param replace.infinite Deprecated, use \code{handle.infinite} instead.
#' @param handle.infinite String specifying how to handle infinite size factors.
#' If \code{"sanitize"}, infinite size factors are replaced with the largest positive factor in \code{size.factors}.
#' This ensures that any normalized values will be, at least, finite; the choice of a relatively large replacement value reflects the extremity of the scaling.
#' If \code{"error"}, an error is raised upon encountering an infinite size factor. 
#' If \code{"ignore"}, no action is taken.
#'
#' If \code{NULL}, the default value in \code{\link{sanitizeSizeFactorsDefaults}} is used.
#' @param replace.nan Deprecated, use \code{handle.nan} instead.
#' @param handle.nan String specifying how to handle NaN size factors.
#' If \code{"sanitize"}, NaN size factors are replaced with unity, e.g., scaling normalization is a no-op.
#' If \code{"error"}, an error is raised upon encountering an NaN size factor. 
#' If \code{"ignore"}, no action is taken.
#' 
#' If \code{NULL}, the default value in \code{\link{sanitizeSizeFactorsDefaults}} is used.
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
sanitizeSizeFactors <- function(
    size.factors,
    replace.zero = NULL,
    replace.negative = NULL,
    replace.infinite = NULL,
    replace.nan = NULL,
    handle.zero = "sanitize",
    handle.negative = "sanitize",
    handle.infinite = "sanitize",
    handle.nan = "sanitize"
) {
    sanitize_size_factors(
        size.factors,
        handle_zero = update_handle(replace.zero, handle.zero, "zero"),
        handle_negative = update_handle(replace.negative, handle.negative, "negative"),
        handle_infinite = update_handle(replace.infinite, handle.infinite, "infinite"),
        handle_nan = update_handle(replace.nan, handle.nan, "nan")
    )
}

update_handle <- function(replace, handle, msg) {
    if (isTRUE(replace)) {
        .Deprecated(old=paste0("replace.", msg, "="), new=paste0("handle.", msg, "="))
        "sanitize"
    } else if (isFALSE(replace)) {
        .Deprecated(old=paste0("replace.", msg, "="), new=paste0("handle.", msg, "="))
        "ignore"
    } else {
        handle
    }
}

#' Default parameters for \code{\link{sanitizeSizeFactors}}
#' @description Default parameters from the underlying C++ library.
#' These may be overridden by defaults in the \code{\link{sanitizeSizeFactors}} function signature.
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' sanitizeSizeFactorsDefaults()
#' @export
sanitizeSizeFactorsDefaults <- function() sanitize_size_factors_defaults()
