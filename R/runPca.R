#' Principal components analysis
#'
#' Run a PCA on the gene-by-cell log-expression matrix and extract the top principal components (PCs).
#' This yields a low-dimensional representation that reduces noise and compute time in downstream analyses.
#' For efficiency, the PCA itself is approximated using IRLBA.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Typically, the matrix is expected to contain log-expression values (see \code{\link{normalizeCounts}}) for \dQuote{interesting} genes (see \code{\link{chooseHighlyVariableGenes}}).
#' @param number Integer scalar specifying the number of top PCs to retain.
#' More PCs will capture more biological signal at the cost of increasing noise and compute time.
#' If this is greater than the maximum number of PCs (i.e., the smaller dimension of \code{x}), only the maximum number of PCs will be reported in the results.
#' @param scale Logical scalar indicating whether to scale all genes to have the same variance.
#' This ensures that each gene contributes equally to the PCA, favoring consistent variation across many genes rather than large variation in a few genes.
#' If \code{block} is specified, each gene's variance is calculated as a weighted sum of the variances from each block. 
#' Genes with zero variance are ignored.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' The PCA will be performed on the residuals after regressing out the block effect, ensuring that differences between block do not dominate the variation in the dataset.
#' Alternatively \code{NULL} if all cells are from the same block.
#' @param block.weight.policy String specifying the policy to use for weighting the contribution of different blocks to the PCA.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL} and \code{block.weight.policy = "variable"}.
#' @param components.from.residuals Logical scalar indicating whether to compute the PC scores from the residuals in the presence of a blocking factor.
#' By default, the residuals are only used to compute the rotation matrix, and the original expression values of the cells are projected onto this new space (see Details).
#' Only used if \code{block} is not \code{NULL}.
#' @param extra.work Integer scalar specifying the extra dimensions for the IRLBA workspace.
#' Larger values improve accuracy at the cost of compute time.
#' @param iterations Integer scalar specifying the maximum number of restart iterations for IRLBA.
#' Larger values improve accuracy at the cost of compute time.
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
#' This can be used to divide \code{variance.explained} to obtain the proportion of variance explained by each PC.
#' \item \code{center}, a numeric vector containing the mean for each gene.
#' If \code{block} is provided, this is instead a matrix containing the mean for each gene (column) in each block (row).
#' \item \code{scale}, a numeric vector containing the scaling for each gene.
#' Only reported if \code{scale=TRUE}.
#' }
#'
#' @details
#' When \code{block} is specified, the nature of the reported PC scores depends on the choice of \code{components.from.residuals}:
#' \itemize{
#' \item If \code{TRUE}, the PC scores are computed from the matrix of residuals.
#' This yields a low-dimensional space where inter-block differences have been removed,
#' assuming that all blocks have the same subpopulation composition and the inter-block differences are consistent for all cell subpopulations.
#' Under these assumptions, we could use these components for downstream analysis without any concern for block-wise effects.
#' \item If \code{FALSE}, the rotation vectors are first computed from the matrix of residuals.
#' To obtain PC scores, each cell is then projected onto the associated subspace using its original expression values.
#' This approach ensures that inter-block differences do not contribute to the PCA but does not attempt to explicitly remove them.
#' }
#' In complex datasets, the assumptions mentioned for \code{TRUE} not hold and more sophisticated batch correction methods like MNN correction are required.
#' Functions like \code{\link{correctMnn}} will accept a low-dimensional embedding of cells that can be created as described above with \code{FALSE}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{simple_pca} and \code{blocked_pca} functions for \url{https://libscran.github.io/scran_pca/}. 
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

    out <- run_pca(
        initializeCpp(x, .check.na=FALSE),
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
        colnames(out$center) <- rownames(x)
    } else {
        names(out$center) <- rownames(x)
    }

    if (!scale) {
        out <- out[setdiff(names(out), "scale")]
    } else {
        names(out$scale) <- rownames(x)
    }

    rownames(out$rotation) <- rownames(x)

    out
}

