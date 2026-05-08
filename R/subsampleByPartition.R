#' Subsample by partition
#'
#' Subsample from a partitioning of cells, preserving the frequency of cells across partitions.
#' The subsample is typically used as a smaller representative dataset that still exhibits the same distribution of cells as the full dataset.
#'
#' @param partitions Factor, integer vector or character vector of length equal to the number of cells,
#' containing the assignment of each cell to a partition (e.g., clustering).
#' @param number Integer specifying the number of cells to retain in the subsample.
#' If greater than the length of \code{partitions}, all cells are retained.
#' @param seed Integer specifying the seed for the random number generator.
#'
#' If \code{NULL}, the default value in \code{\link{subsampleByartitionDefaults}} is used.
#' @param force.non.empty Boolean indicating whether each partition should have at least one cell in the subsample.
#' If \code{FALSE}, partitions may not be represented if the number of cells is less than the ratio of \code{number} to \code{length(partitions)}.
#'
#' If \code{NULL}, the default value in \code{\link{subsampleByartitionDefaults}} is used.
#'
#' @return Integer vector containing the indices of the cells retained in the subsample.
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{subsampleByNeighbors}}, for an alternative downsampling strategy.
#' 
#' @examples
#' x <- matrix(rnorm(10000), nrow=5)
#' partitions <- clusterKmeans(x, k=20)
#' keep <- subsampleByPartition(partitions$clusters, 100)
#' table(partitions$clusters[keep])
#' sub.x <- x[,keep] # subsample to be used for expensive downstream steps.
#'
#' # Rare partitions are retained by default:
#' grouping <- c(LETTERS[1:10], sample(LETTERS[11:26], 990, replace=TRUE))
#' keep <- subsampleByPartition(grouping, 100)
#' table(grouping[keep])
#' keep <- subsampleByPartition(grouping, 100, force.non.empty=FALSE)
#' table(grouping[keep])
#'
#' @export
subsampleByPartition <- function(partitions, number, seed = NULL, force.non.empty = NULL) { 
    part <- .transformFactor(partitions)
    subsample_by_partition(part$index, number, seed = seed, force_non_empty = force.non.empty)
}

#' Default parameters for \code{\link{subsampleByPartition}}
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' subsampleByPartitionDefaults()
#' @export
subsampleByPartitionDefaults <- function() subsample_by_partition_defaults()
