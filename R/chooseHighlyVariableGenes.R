#' Choose highly variable genes
#'
#' Choose highly variable genes based on a variance-related statistic.
#'
#' @param stats Numeric vector of variances (or a related statistic) across all genes.
#' Typically the residuals from \code{\link{modelGeneVariances}} are used here.
#' @param top Integer specifying the number of top genes to retain.
#' @param larger Logical scalar indicating whether larger values of \code{stats} correspond to more variable genes.
#' @param keep.ties Logical scalar indicating whether to keep tied values of \code{stats}, even if \code{top} may be exceeded.
#'
#' @return Integer vector containing the indices of genes in \code{stats} that are considered to be highly variable.
#' 
#' @examples
#' resids <- rexp(10000)
#' str(chooseHighlyVariableGenes(resids))
#'
#' @seealso
#' \url{https://libscran.github.io/scran_variances/}, for the underlying implementation.
#'
#' @author Aaron Lun
#' @export
chooseHighlyVariableGenes <- function(stats, top=4000, larger=TRUE, keep.ties=TRUE) {
    out <- choose_highly_variable_genes(stats, top=top, larger=larger, keep_ties=keep.ties)
    out + 1L
}
