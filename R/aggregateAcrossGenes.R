#' Aggregate expression across genes
#'
#' Aggregate expression values across genes, potentially with weights.
#' This is typically used to summarize expression values for gene sets into a single per-cell score.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Values are usually normalized expression values, possibly log-transformed depending on the application.
#' @param sets A list of integer vectors containing the row indices of genes in each set.
#' Alternatively, each entry may be a list of length 2, containing an integer vector (row indices) and a numeric vector (weights).
#' @param average Logical scalar indicating whether to compute the average rather than the sum.
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list of length equal to that of \code{sets}.
#' Each entry is a numeric vector of length equal to the number of columns in \code{x}, 
#' containing the (weighted) sum/mean of expression values for the corresponding set across all cells.
#'
#' @author Aaron Lun
#' @seealso
#' The \code{aggregate_across_genes} function in \url{https://libscran.github.io/scran_aggregate/}. 
#'
#' \code{\link{aggregateAcrossCells}}, to aggregate expression values across groups of cells.
#'
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Unweighted aggregation:
#' sets <- list(
#'    foo = sample(nrow(x), 20),
#'    bar = sample(nrow(x), 10)
#' )
#' agg <- aggregateAcrossGenes(x, sets)
#' str(agg)
#'
#' # Weighted aggregation:
#' sets <- list(
#'    foo = list(sample(nrow(x), 20), runif(20)),
#'    bar = list(sample(nrow(x), 10), runif(10))
#' )
#' agg2 <- aggregateAcrossGenes(x, sets, average = TRUE)
#' str(agg2)
#' 
#' @export
#' @importFrom beachmat initializeCpp
aggregateAcrossGenes <- function(x, sets, average = FALSE, num.threads = 1) {
    ptr <- initializeCpp(x, .check.na=FALSE)
    output <- aggregate_across_genes(ptr, sets, average, num.threads)
    names(output) <- names(sets)
    output
}
