#' Score gene set activity for each cell
#'
#' Compute per-cell scores for a gene set, defined as the column sums of a rank-1 approximation to the submatrix for the gene set.
#' This uses the same approach as the \pkg{GSDecon} package by Jason Hackney, adapted to use an approximate PCA (via IRLBA) and to support blocking.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Typically, the matrix is expected to contain log-expression values.
#' @param set Vector specifying the rows of \code{x} that belong to the gene set.
#' This may be an integer vector of row indices, a logical vector of length equal to the number of rows, or a character vector of row names.
#' For integer and character vectors, duplicate entries are ignored.
#' @param rank Integer scalar specifying the rank of the approximation.
#' The default value of 1 assumes that each gene set only describes a single coordinated biological function.
#' @inheritParams runPca
#'
#' @return List containing:
#' \itemize{
#' \item \code{scores}, a numeric vector of per-cell scores for each column in \code{x}.
#' \item \code{weights}, a data frame containing \code{row}, an integer vector of ordered and unique row indices corresponding to the genes in \code{set};
#' and \code{weight}, a numeric vector of per-gene weights for each gene in \code{row}. 
#' }
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{compute} and \code{compute_blocked} functions in \url{https://libscran.github.io/gsdecon/}.
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
    chosen <- .sanitizeGeneSet(set, n=nr, names=rownames(x))
    ptr <- tatami.subset(ptr, chosen, by.row=TRUE)

    out <- score_gene_set(
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

    out$weights <- data.frame(row=chosen, weight=out$weights)
    out
}

.sanitizeGeneSet <- function(set, n, names, arg = "set") {
    # Make sure we obtain unique and sorted indices so that
    # the slicing operations in tatami.subset are efficient.
    if (is.numeric(set)) {
        set <- as.integer(set)
        if (anyDuplicated(set)) {
            set <- unique(set)
        }
        if (anyNA(set) || min(set) < 1 || max(set) > n) {
            stop("'", arg, "' contains out-of-range indices")
        }
        if (is.unsorted(set)) {
            set <- sort(set)
        }
        return(set)

    } else if (is.logical(set)) {
        if (length(set) != n) {
            stop("length of '", arg, "' should be equal to the number of rows")
        }
        return(which(set))

    } else {
        if (anyDuplicated(set)) {
            set <- unique(set)
        }
        set <- match(set, names)
        if (anyNA(set)) {
            set <- set[!is.na(set)]
        }
        if (is.unsorted(set)) {
            set <- sort(set)
        }
        return(set)
    }
}
