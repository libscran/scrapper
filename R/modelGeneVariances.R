#' Model per-gene variances in expression
#'
#' Compute the variance in (log-)expression values for each gene,
#' and model the trend in the variances with respect to the mean.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' It is typically expected to contain log-expression values, e.g., from \code{\link{normalizeCounts}}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' Alternatively \code{NULL} if all cells are from the same block.
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic
#' (See documentation at \url{https://libscran.github.io/scran_blocks} for more details.)
#' Only used if \code{block} is not \code{NULL}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' The first and second values are used as the lower and upper bounds, respectively, for the variable weight calculation.
#' (See documentation at \url{https://libscran.github.io/scran_blocks} for more details.)
#' Only used if \code{block} is not \code{NULL} and \code{block.weight.policy = "variable"}.
#' @inheritParams fitVarianceTrend
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return A list containing \code{means}, \code{variances}, \code{fitted} and \code{residuals},
#' each of which is a numeric vector containing the statistic of the same name across all genes.
#'
#' If \code{block} is supplied, each of the vectors contains the average across all blocks.
#' The list will also contain \code{per.block}, a list of lists containing the equivalent statistics for each block.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{model_gene_variances} function in \url{https://libscran.github.io/scran_variances}.
#'
#' \code{\link{fitVarianceTrend}}, which fits the mean-variance trend.
#'
#' @examples
#' library(Matrix)
#' x <- abs(rsparsematrix(1000, 100, 0.1) * 10)
#' out <- modelGeneVariances(x)
#' str(out)
#'
#' # Throwing in some blocking.
#' block <- sample(letters[1:4], ncol(x), replace=TRUE)
#' out <- modelGeneVariances(x, block=block)
#' str(out)
#'
#' @export
#' @importFrom beachmat initializeCpp
modelGeneVariances <- function(
    x,
    block=NULL,
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    mean.filter=TRUE,
    min.mean=0.1, 
    transform=TRUE, 
    span=0.3,
    use.min.width=FALSE,
    min.width=1,
    min.window.count=200,
    num.threads=1) 
{
    block <- .transformFactor(block)

    output <- model_gene_variances(
        initializeCpp(x),
        block=block$index,
        nblocks=length(block$names),
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        mean_filter=mean.filter,
        min_mean=min.mean,
        transform=transform,
        span=span,
        use_min_width=use.min.width,
        min_width=min.width,
        min_window_count=min.window.count,
        num_threads=num.threads
    )

    if (!is.null(block$index)) {
        names(output$per.block) <- block$names
    }

    output
}
