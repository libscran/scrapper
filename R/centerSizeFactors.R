#' Center size factors
#'
#' Scale the size factors so they are centered at unity,
#' which ensures that the scale of the counts is preserved (on average) after normalization.
#'
#' @param size.factors Numeric vector of size factors across cells.
#' @param block Vector or factor of length equal to \code{size.factors}, specifying the block of origin for each cell.
#' Alternatively \code{NULL}, in which case all cells are assumed to be in the same block.
#' @param mode String specifying how to scale size factors across blocks.
#' \code{"lowest"} will compute the average size factor in each block, identify the lowest average across all blocks, and then scale all size factors by that value.
#' \code{"per-block"} will compute the average size factor in each block, and then scale each size factor by the average of block to which it belongs.
#' Only used if \code{block} is provided.
#' 
#' @return Numeric vector of length equal to \code{size.factors}, containing the centered size factors.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/scran_norm/}, for the rationale behind centering the size factors.
#'
#' @examples
#' centerSizeFactors(runif(100))
#'
#' centerSizeFactors(runif(100), block=sample(3, 100, replace=TRUE))
#'
#' @export
centerSizeFactors <- function(size.factors, block=NULL, mode=c("lowest", "per-block")) {
    block <- .transformFactor(block)
    center_size_factors(size.factors, block$index, match.arg(mode) == "lowest")
}
