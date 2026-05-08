#' Choose highly variable genes
#'
#' Choose highly variable genes (HVGs) based on a variance-related statistic.
#' This is typically used to subset the gene-cell matrix prior to calling \code{\link{runPca}}.
#'
#' @param stats Numeric vector of variances (or a related statistic) across all genes.
#' Typically, the residuals from \code{\link{modelGeneVariances}} are used here.
#' @param top Integer specifying the number of top genes to retain.
#' Note that the actual number of retained genes may not be equal to \code{top}, depending on the other options.
#'
#' If \code{NULL}, the default value in \code{\link{chooseHighlyVariableGenesDefaults}} is used.
#' @param larger Boolean indicating whether larger values of \code{stats} correspond to more variable genes.
#' If \code{TRUE}, HVGs are defined as those with the largest values of \code{stats}.
#' This is typically the case for variances or related statistics, e.g., residuals.
#'
#' If \code{NULL}, the default value in \code{\link{chooseHighlyVariableGenesDefaults}} is used.
#' @param keep.ties Boolean indicating whether to keep tied values of \code{stats}, even if \code{top} may be exceeded.
#'
#' If \code{NULL}, the default value in \code{\link{chooseHighlyVariableGenesDefaults}} is used.
#' @param bound Number specifying the lower bound (if \code{larger=TRUE}) or upper bound (otherwise) to be applied to \code{stats}.
#' Genes are not considered to be HVGs if they do not satisfy this bound, even if they are within the \code{top} genes.
#' For example, residuals from the fitted trend should be positive, which can be enforced by setting \code{bound = 0}. 
#'
#' If \code{NULL}, the default value in \code{\link{chooseHighlyVariableGenesDefaults}} is used.
#'
#' If \code{NA}, no bound is enforced.
#'
#' @return Integer vector containing the indices of genes in \code{stats} that are considered to be highly variable.
#' 
#' @examples
#' resids <- rexp(10000)
#' str(chooseHighlyVariableGenes(resids))
#'
#' @seealso
#' The \code{choose_highly_variable_genes} function in \url{https://libscran.github.io/scran_variances/}. 
#'
#' \code{\link{chooseRnaHvgs.se}}, which choose the HVGs from the residuals computed from a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @author Aaron Lun
#' @export
chooseHighlyVariableGenes <- function(stats, top = NULL, larger = NULL, keep.ties = NULL, bound = NULL) {
    if (!is.null(top)) {
        top <- min(length(stats), top) # protect against top=Inf, which Rcpp can't cast to an int.
    }
    out <- choose_highly_variable_genes(stats, top=top, larger=larger, keep_ties=keep.ties, bound=bound)
    out + 1L
}

#' Default parameters for \code{\link{chooseHighlyVariableGenes}}
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' chooseHighlyVariableGenesDefaults()
#' @export
chooseHighlyVariableGenesDefaults <- function() choose_highly_variable_genes_defaults()
