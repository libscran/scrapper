#' Aggregate expression across cells
#'
#' Aggregate expression values across cells based on one or more grouping factors.
#' This is primarily used to create pseudo-bulk profiles for each cluster/sample combination.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Values are typically expected to be counts.
#' @param factors A list or data frame containing one or more grouping factors, see \code{\link{combineFactors}}.
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sums}, a numeric matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the summed expression across all cells with that combination. 
#' \item \code{detected}, an integer matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the number of cells with detected expression in that combination.
#' \item \code{combinations}, a data frame describing the levels for each unique combination of factors.
#' Rows of this data frame correspond to columns of \code{sums} and \code{detected},
#' while columns correspond to the factors in \code{factors}.
#' \item \code{counts}, the number of cells associated with each combination.
#' Each entry corresponds to a row of \code{combinations}.
#' \item \code{index}, an integer vector of length equal to the number of cells in \code{x}.
#' This specifies the combination in \code{combinations} to which each cell was assigned.
#' }
#'
#' @seealso
#' \code{\link{aggregateAcrossGenes}}, to aggregate expression values across gene sets.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Simple aggregation:
#' clusters <- sample(LETTERS, 100, replace=TRUE)
#' agg <- aggregateAcrossCells(x, list(cluster=clusters))
#' str(agg)
#'
#' # Multi-factor aggregation
#' samples <- sample(1:5, 100, replace=TRUE)
#' agg2 <- aggregateAcrossCells(x, list(cluster=clusters, sample=samples))
#' str(agg2)
#' 
#' @export
#' @importFrom beachmat initializeCpp
aggregateAcrossCells <- function(x, factors, num.threads = 1) {
    combined <- combineFactors(factors)

    ptr <- initializeCpp(x)
    output <- aggregate_across_cells(ptr, combined$index - 1L, num.threads)
    rownames(output$sums) <- rownames(output$detected) <- rownames(x)

    output$combinations <- combined$levels
    output$counts <- tabulate(combined$index, nbins=nrow(output$combinations))
    output$index <- combined$index

    output
}
