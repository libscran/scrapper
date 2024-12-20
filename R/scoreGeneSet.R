#' Score gene set activity for each cell
#'
#' Compute per-cell scores for a gene set, defined as the column sums of a rank-1 approximation to the submatrix for the feature set.
#' This uses the same approach implemented in the \pkg{GSDecon} package from Jason Hackney.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Typically, the matrix is expected to contain log-expression values, and the rows should be filtered to relevant (e.g., highly variable) genes.
#' @param set Integer, logical or character vector specifying the rows that belong to the gene set.
#' @param rank Integer scalar specifying the rank of the approximation.
#' @inheritParams runPca
#'
#' @return List containing \code{scores}, a numeric vector of per-cell scores for each column in \code{x};
#' and \code{weights}, a numeric vector of per-feature weights for each feature in \code{set}. 
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/gsdecon/}, for more details on the underlying algorithm.
#'
#' @author Aaron Lun
#' @examples
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' normed <- normalizeCounts(x, size.factors=centerSizeFactors(colSums(x)))
#' scoreGeneSet(normed, set=c(1,3,5,10,20,100))
#' 
#' @export
#' @importFrom beachmat initializeCpp tatami.dim tatami.subset
scoreGeneSet <- function(
    x,
    set,
    rank=1, 
    scale=FALSE,
    block=NULL, 
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    extra.work=7,
    iterations=1000,
    seed=5489,
    realized=TRUE,
    num.threads=1) 
{
    block <- .transformFactor(block)
    ptr <- initializeCpp(x, .check.na=FALSE)

    nr <- tatami.dim(ptr)[1]
    chosen <- which(.toLogical(set, n=nr, names=rownames(x)))
    ptr <- tatami.subset(ptr, chosen, by.row=TRUE)

    score_gene_set(
        ptr,
        block=block$index,
        rank=rank,
        scale=scale,
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        realized=realized,
        irlba_work=extra.work,
        irlba_iterations=iterations,
        irlba_seed=seed,
        num_threads=num.threads
    )
}
