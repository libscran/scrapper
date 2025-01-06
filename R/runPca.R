#' Principal components analysis
#'
#' Run a PCA on the gene-by-cell log-expression matrix to obtain a low-dimensional representation for downstream analyses.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Typically, the matrix is expected to contain log-expression values, and the rows should be filtered to relevant (e.g., highly variable) genes.
#' @param number Integer scalar specifying the number of PCs to retain.
#' @param scale Logical scalar indicating whether to scale all genes to have the same variance.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' Alternatively \code{NULL} if all cells are from the same block.
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic.
#' Only used if \code{block} is not \code{NULL}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' The first and second values are used as the lower and upper bounds, respectively, for the variable weight calculation.
#' Only used if \code{block} is not \code{NULL} and \code{block.weight.policy = "variable"}.
#' @param components.from.residuals Logical scalar indicating whether to compute the PC scores from the residuals in the presence of a blocking factor.
#' By default, the residuals are only used to compute the rotation matrix, and the original expression values of the cells are projected onto this new space.
#' Only used if \code{block} is not \code{NULL}.
#' @param extra.work Integer scalar specifying the extra dimensions for the IRLBA workspace.
#' @param iterations Integer scalar specifying the maximum number of restart iterations for IRLBA.
#' @param seed Integer scalar specifying the seed for the initial random vector in IRLBA.
#' @param realized Logical scalar indicating whether to realize \code{x} into an optimal memory layout for IRLBA.
#' This speeds up computation at the cost of increased memory usage.
#' @param num.threads Number of threads to use.
#'
#' @return List containing:
#' \itemize{
#' \item \code{components}, a matrix of PC scores.
#' Rows are dimensions (i.e., PCs) and columns are cells.
#' \item \code{rotation}, the rotation matrix.
#' Rows are genes and columns are dimensions.
#' \item \code{variance.explained}, the vector of variances explained by each PC.
#' \item \code{total.variance}, the total variance in the dataset.
#' \item \code{center}, a numeric vector containing the mean for each gene.
#' If \code{block} is provided, this is instead a matrix containing the mean for each gene (column) in each block (row).
#' \item \code{scale}, a numeric vector containing the scaling for each gene.
#' Only reported if \code{scale=TRUE}.
#' }
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/scran_pca/}, for more details on the PCA.
#' In particular, the documentation for the \code{blocked_pca} function explains the blocking strategy.
#'
#' @examples
#' library(Matrix)
#' x <- abs(rsparsematrix(1000, 100, 0.1) * 10)
#' y <- normalizeCounts(x, size.factors=centerSizeFactors(colSums(x)))
#'
#' # A simple PCA:
#' out <- runPca(y)
#' str(out)
#'
#' # Blocking on uninteresting factors:
#' block <- sample(LETTERS[1:3], ncol(y), replace=TRUE)
#' bout <- runPca(y, block=block)
#' str(bout)
#' 
#' @export
#' @importFrom beachmat initializeCpp
runPca <- function(x,
    number=25,
    scale=FALSE,
    block=NULL, 
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    components.from.residuals=FALSE,
    extra.work=7,
    iterations=1000,
    seed=5489,
    realized=TRUE,
    num.threads=1) 
{
    block <- .transformFactor(block)

    x <- initializeCpp(x, .check.na=FALSE)
    out <- run_pca(x, 
        number=number,
        scale=scale,
        block=block$index,
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        components_from_residuals=components.from.residuals,
        realized=realized,
        irlba_work=extra.work,
        irlba_iterations=iterations,
        irlba_seed=seed,
        num_threads=num.threads
    )

    if (!is.null(block$index)) {
        rownames(out$center) <- block$names
    }

    if (!scale) {
        out <- out[setdiff(names(out), "scale")]
    }

    out
}

