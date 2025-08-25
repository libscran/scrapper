#' Choose a suitable pseudo-count
#'
#' Choose a suitable pseudo-count to control the bias introduced by log-transformation of normalized counts from \code{\link{normalizeCounts}}.
#' Larger pseudo-counts shrink log-expression values towards the zero-expression baseline,
#' reducing the impact of the transformation bias at the cost of some sensitivity.
#'
#' @param size.factors Numeric vector of size factors for all cells.
#' It is expected that these have already been centered, e.g., with \code{\link{centerSizeFactors}}.
#' Invalid size factors (e.g., non-positive, non-finite) will be ignored.
#' @param quantile Numeric scalar specifying the quantile to use for finding the smallest/largest size factors.
#' Setting this to zero will use the observed minimum and maximum, though in practice, this is usually too sensitive to outliers.
#' The default is to take the 5th and 95th percentile to obtain a range that captures most of the distribution.
#' @param max.bias Numeric scalar specifying the maximum allowed bias.
#' This is the maximum absolute value of any spurious log2-fold change between the cells with the smallest and largest size factors.
#' @param min.value Numeric scalar specifying the minimum value for the pseudo-count.
#' Defaults to 1 to stabilize near-zero normalized expression values, otherwise these manifest as avoid large negative values.
#'
#' @return A choice of pseudo-count for \code{\link{normalizeCounts}}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{choose_pseudo_count} in \url{https://libscran.github.io/scran_norm/}.
#'
#' Lun ATL (2018).
#' Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data.
#' \emph{biorXiv} doi:10.1101/404962
#'
#' @examples
#' sf <- centerSizeFactors(runif(100))
#' choosePseudoCount(sf)
#' choosePseudoCount(sf, quantile=0.01)
#' choosePseudoCount(sf, max.bias=0.5)
#'
#' @export
choosePseudoCount <- function(size.factors, quantile = 0.05, max.bias = 1, min.value = 1) {
    choose_pseudo_count(size.factors, quantile, max.bias, min.value)
}
