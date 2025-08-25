#' Model per-gene variances in expression
#'
#' Model the per-gene variances as a function of the mean in single-cell expression data.
#' Highly variable genes can then be selected for downstream analyses. 
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' It is typically expected to contain log-expression values, e.g., from \code{\link{normalizeCounts}}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' If provided, calculation of means/variances and trend fitting are performed within each block to ensure that block effects do not confound the estimates.
#' Alternatively, \code{NULL} if all cells are from the same block.
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic.
#' This should be one of:
#' \itemize{
#' \item \code{"none"}: the contribution of each block is proportional to its size.
#' \item \code{"equal"}: blocks are equally weighted regardless of their size.
#' \item \code{"variable"}: blocks are equally weighted past a certain threshold size.
#' Below that size, the contribution of each block is proportional to its size.
#' This avoids outsized contributions from very large blocks.
#' }
#' Only used if \code{block} is not \code{NULL}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' The first value is usually zero and defines the threshold on the size at or below which a block receives zero weight.
#' The second value is the upper threshold on the size above which all blocks have the same weight.
#' Only used if \code{block} is not \code{NULL} and \code{block.weight.policy = "variable"}.
#' @inheritParams fitVarianceTrend
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @details
#' We compute the mean and variance for each gene and fit a trend to the variances with respect to the means using \code{\link{fitVarianceTrend}}.
#' We assume that most genes at any given abundance are not highly variable, such that the fitted value of the trend is interpreted as the \dQuote{uninteresting} variance - 
#' this is mostly attributed to technical variation like sequencing noise, but can also represent constitutive biological noise like transcriptional bursting.
#' Under this assumption, the residual can be treated as a measure of biologically interesting variation.
#' Genes with large residuals can then be selected for downstream analyses, e.g., with \code{\link{chooseHighlyVariableGenes}}.
#'
#' @return A list containing \code{statistics}, a data frame with number of rows equal to the number of genes.
#' This contains the columns \code{means}, \code{variances}, \code{fitted} and \code{residuals},
#' each of which is a numeric vector containing the statistic of the same name across all genes.
#'
#' If \code{block} is supplied, each of the column vectors described above contains the average across all blocks.
#' The list will also contain \code{per.block}, a list of data frames containing the equivalent statistics for each block.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{model_gene_variances} function in \url{https://libscran.github.io/scran_variances/}.
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

    stats <- model_gene_variances(
        initializeCpp(x, .check.na=FALSE),
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

    output <- list(
        statistics = data.frame(
            means = stats$means,
            variances = stats$variances,
            fitted = stats$fitted,
            residuals = stats$residuals,
            row.names = rownames(x)
        )
    )

    if (!is.null(block$index)) {
        pb <- stats$per.block
        for (i in seq_along(pb)) {
            curstats <- pb[[i]]
            pb[[i]] <- data.frame(
                means = curstats$means,
                variances = curstats$variances,
                fitted = curstats$fitted,
                residuals = curstats$residuals,
                row.names = rownames(x)
            )
        }
        names(pb) <- block$names
        output$per.block <- pb
    }

    output
}
