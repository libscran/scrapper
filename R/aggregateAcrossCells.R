#' Aggregate expression across cells
#'
#' Aggregate expression values across cells based on one or more grouping factors.
#' This is usually applied to a count matrix to create pseudo-bulk profiles for each cluster/sample combination.
#' These profiles can then be used as if they were counts from bulk data, e.g., for differential analyses with \pkg{edgeR}.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Values are typically expected to be counts.
#' @param factors A list or data frame (or their equivalents from \pkg{S4Vectors}) containing one or more grouping factors, see \code{\link{combineFactors}}.
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sums}, a numeric matrix where each row corresponds to a gene and each column corresponds to a unique combination of levels from \code{factors}.
#' Each entry contains the summed expression across all cells with that combination. 
#' \item \code{detected}, an integer matrix where each row corresponds to a gene and each column corresponds to a unique combination of levels from \code{factors}.
#' Each entry contains the number of cells with detected expression in that combination.
#' \item \code{combinations}, a \link[S4Vectors]{DataFrame} containing the unique combination of levels from \code{factors}.
#' Rows correspond to columns of \code{sums} and \code{detected}, while columns correspond to the factors in \code{factors}.
#' \item \code{counts}, the number of cells associated with each combination.
#' Each entry corresponds to a row of \code{combinations}.
#' \item \code{index}, an integer vector of length equal to the number of cells in \code{x}.
#' This specifies the combination in \code{combinations} to which each cell was assigned.
#' }
#'
#' @seealso
#' The \code{aggregate_across_cells} function in \url{https://libscran.github.io/scran_aggregate/}.
#'
#' \code{\link{aggregateAcrossCells.se}}, to perform aggregation on a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
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
#' @importFrom S4Vectors DataFrame
aggregateAcrossCells <- function(x, factors, num.threads = 1) {
    .checkSEX(x, "aggregateAcrossCells.se")

    combined <- combineFactors(factors)

    ptr <- initializeCpp(x, .check.na=FALSE)
    output <- aggregate_across_cells(ptr, combined$index - 1L, num.threads)
    rownames(output$sums) <- rownames(output$detected) <- rownames(x)

    output$combinations <- combined$levels
    output$counts <- tabulate(combined$index, nbins=nrow(output$combinations))
    output$index <- combined$index

    output
}
