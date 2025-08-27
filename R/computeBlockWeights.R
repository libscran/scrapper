#' Compute block weights
#'
#' Compute a weight for each block based on the number of cells in each block.
#' This is typically used to aggregate statistics across blocks, e.g., with weighted sums/averages.
#' 
#' @param sizes Numeric vector containing the size of (i.e., number of cells in) each block.
#' @param block.weight.policy String specifying the policy to use for weighting different blocks. 
#' This should be one of:
#' \itemize{
#' \item \code{"size"}: the contribution of each block is proportional to its size.
#' \code{"none"} is also a deprecated alias for \code{"size"}. 
#' \item \code{"equal"}: blocks are equally weighted regardless of their size.
#' The exception is that of empty blocks with no cells, which receive zero weight.
#' \item \code{"variable"}: blocks are equally weighted past a certain threshold size.
#' Below that size, the contribution of each block is proportional to its size.
#' This avoids outsized contributions from very large blocks.
#' }
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' The first and second values are used as the lower and upper bounds, respectively, for the variable weight calculation.
#' Only used if \code{block.weight.policy = "variable"}.
#'
#' @return Numeric vector containing the relative block weights.
#'
#' @author Aaron Lun
#' @examples
#' computeBlockWeights(c(1, 10, 100, 1000, 10000))
#' computeBlockWeights(c(1, 10, 100, 1000, 10000), block.weight.policy="equal")
#' computeBlockWeights(c(1, 10, 100, 1000, 10000), variable.block.weight=c(50, 5000))
#'
#' @seealso
#' The \code{compute_weights} function from \url{https://libscran.github.io/scran_blocks/}.
#'
#' @export
computeBlockWeights <- function(
    sizes,
    block.weight.policy=c("variable", "equal", "size", "none"),
    variable.block.weight=c(0, 1000))
{
    compute_block_weights(sizes, match.arg(block.weight.policy), variable.block.weight)
}
