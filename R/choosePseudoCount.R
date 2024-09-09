#' Choose a suitable pseudo-count
#'
#' Choose a suitable pseudo-count to control the bias introduced by log-transformation of normalized counts.
#'
#' @param size.factors Numeric vector of size factors for all cells.
#' @param quantile Numeric scalar specifying the quantile to use for defining extreme size factors.
#' @param max.bias Numeric scalar specifying the maximum allowed bias.
#' @param min.value Numeric scalar specifying the minimum value for the pseudo-count.
#'
#' @return A choice of pseudo-count for \code{\link{normalizeCounts}}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/scran_norm/}, for the motivation behind calculating a larger pseudo-count.
#'
#' @examples
#' sf <- runif(100)
#' choosePseudoCount(sf)
#' choosePseudoCount(sf, quantile=0.01)
#' choosePseudoCount(sf, max.bias=0.5)
#'
#' @export
choosePseudoCount <- function(size.factors, quantile = 0.05, max.bias = 1, min.value = 1) {
    choose_pseudo_count(size.factors, quantile, max.bias, min.value)
}
