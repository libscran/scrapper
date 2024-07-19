#' Combine multiple factors 
#'
#' Combine multiple categorical factors based on the unique combinations of levels from each factor.
#'
#' @param factors List of vectors or factors of the same length.
#' Corresponding elements across all vectors/factors represent the combination of levels for a single observation.
#' For factors, any existing levels are respected.
#' For other vectors, the sorted and unique values are used as levels.
#'
#' Alternatively, a data frame where each column is a vector or factor and each row corresponds to an observation.
#' @param keep.unused Logical scalar indicating whether to report unused combinations of levels.
#'
#' @return List containing \code{levels}, a data frame containing the sorted and unique combinations of levels from \code{factors};
#' and \code{index}, an integer vector specifying the index into \code{levels} for each observation.
#'
#' @author Aaron Lun
#'
#' @examples
#' combineFactors(list(
#'     sample(LETTERS[1:5], 100, replace=TRUE),
#'     sample(3, 100, replace=TRUE)
#' ))
#'
#' combineFactors(list(
#'     factor(sample(LETTERS[1:5], 10, replace=TRUE), LETTERS[1:5]),
#'     factor(sample(5, 10, replace=TRUE), 1:5)
#' ), keep.unused=TRUE)
#'
#' @export
combineFactors <- function(factors, keep.unused=FALSE) {
    out <- .combineFactors(factors, keep.unused=keep.unused)
    out$index <- out$index + 1L
    out
}

.combineFactors <- function(factors, keep.unused=FALSE) {
    nfac <- length(factors)
    f0 <- vector("list", nfac)
    levels0 <- vector("list", nfac)

    for (f in seq_len(nfac)) {
        current <- factors[[f]]
        if (is.factor(current)) {
            f0[[f]] <- as.integer(current)
            levels0[[f]] <- levels(current)
        } else {
            curlevels <- sort(unique(current))
            levels0[[f]] <- curlevels
            f0[[f]] <- match(current, curlevels)
        }

        f0[[f]] <- f0[[f]] - 1L # get to 0-based indices.
    }

    combined <- combine_factors(f0, keep.unused, lengths(levels0))

    for (f in seq_len(nfac)) {
        combined$levels[[f]] <- levels0[[f]][combined$levels[[f]] + 1L]
    }
    names(combined$levels) <- names(factors)

    if (is.null(names(combined$levels))) {
        names(combined$levels) <- make.names(seq_len(nfac))
    }
    combined$levels <- data.frame(combined$levels)

    combined 
}
