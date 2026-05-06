#' Center spike-in size factors
#'
#' Center size factors for endogenous genes and spike-in transcripts, following the rationale in \code{\link{centerSizeFactors}}.
#' This aims to preserve the scale of the counts in each feature set so that the normalized expression values are still comparable.
#' 
#' @param endogenous Numeric vector of size factors for endogenous genes across cells.
#' Invalid size factors (e.g., non-positive, non-finite) will be ignored.
#' @param spike.ins List of numeric vectors.
#' Each vector should be of length equal to \code{endogenous} and contains size factors for a set of spike-in transcripts.
#' Invalid size factors (e.g., non-positive, non-finite) will be ignored.
#' @param block Vector or factor of length equal to \code{endogenous}, specifying the block of origin for each cell.
#' Alternatively \code{NULL}, in which case all cells are assumed to be in the same block.
#' @param mode String specifying how to scale size factors across blocks, see the argument of the same name in \code{\link{centerSizeFactors}}.
#'
#' If \code{NULL}, the default value in \code{\link{centerSpikeInFactors}} is used.
#'
#' This argument is only used if \code{block} is provided.
#' 
#' @return List containing:
#' \itemize{
#' \item \code{endogenous}, a numeric vector containing the centered size factors for the endogenous genes.
#' \item \code{spike.ins}, a named list of numeric vectors.
#' Each vector is named after an entry of \code{spike.ins} and contains centered size factors for the corresponding spike-in transcripts.
#' }
#'
#' @details
#' This function is effectively a convenient wrapper around \code{\link{centerSizeFactors}},
#' configured to ensure that the mean of each set of spike-in size factors is the same as that of the endogenous size factors.
#' The aim is to preserve the scale of the average abundances of the spike-ins relative to the endogenous genes for a valid comparison between feature sets.
#' Typically, we would use spike-in transcripts to infer some properties of endogenous genes of similar molarity, e.g., technical variance.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{center_spike_in_factors} and \code{center_spike_in_factors_blocked} functions in \url{https://libscran.github.io/scran_norm/}.
#'
#' @examples
#' endogenous <- runif(100)
#' spike.ercc <- runif(100)
#' centerSpikeInFactors(endogenous, list(ERCC = spike.ercc))
#'
#' block <- sample(3, 100, replace=TRUE)
#' centerSpikeInFactors(endogenous, list(ERCC = spike.ercc), block=block)
#'
#' @export
centerSpikeInFactors <- function(endogenous, spike.ins, block = NULL, mode = NULL) {
    block <- .transformFactor(block)
    output <- center_spike_in_factors(endogenous, spike.ins, block$index, mode)
    names(output$spike.ins) <- names(spike.ins)
    output
}

#' Default parameters for \code{\link{centerSpikeInFactors}}
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' centerSpikeInFactorsDefaults()
#' @export
centerSpikeInFactorsDefaults <- function() center_spike_in_factors_defaults()
