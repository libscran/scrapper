#' Defunct functions
#'
#' See Details for each function's recommended replacements.
#'
#' @param ... Arguments, ignored.
#'
#' @details
#' \code{analyze} is replaced by \code{\link{analyze.se}}.
#'
#' \code{convertAnalyzeResults} is no longer necessary as conversion is done within \code{\link{analyze.se}}.
#'
#' @export
#' @rdname defunct
analyze <- function(...) {
    .Defunct(new="analyze.se")
}

#' @export
#' @rdname defunct
convertAnalyzeResults <- function(...) {
    .Defunct(new="analyze.se")
}
