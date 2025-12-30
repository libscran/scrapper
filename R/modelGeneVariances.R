#' Model per-gene variances in expression
#'
#' Model the per-gene variances as a function of the mean in single-cell expression data.
#' Highly variable genes can then be selected for downstream analyses. 
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' It is typically expected to contain log-expression values, e.g., from \code{\link{normalizeCounts}}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' If provided, calculation of means/variances and trend fitting are performed within each block to ensure that block effects do not confound the estimates.
#' The weighted average of each statistic across all blocks is reported for each gene.
#' Alternatively \code{NULL}, if all cells are from the same block.
#' @param block.average.policy String specifying the policy to use for average statistics across blocks.
#' This can either be a (weighted) \code{"mean"} or a \code{"quantile"}.
#' Only used if \code{block} is not \code{NULL}. 
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL} and \code{block.average.policy="mean"}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL}, \code{block.average.policy="mean"} and \code{block.weight.policy = "variable"}.
#' @param block.quantile Number specifying the probability of the quantile of statistics across blocks. 
#' Defaults to 0.5, i.e., the median of per-block statistics.
#' Only used if \code{block} is not \code{NULL} and \code{block.average.policy="quantile"}.
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
#' @return A list containing \code{statistics}, a \link[S4Vectors]{DataFrame} with number of rows equal to the number of genes.
#' This contains the columns \code{means}, \code{variances}, \code{fitted} and \code{residuals},
#' each of which is a numeric vector containing the statistic of the same name across all genes.
#'
#' If \code{block} is supplied, each of the column vectors described above contains the average across all blocks.
#' The list will also contain \code{per.block}, a list of DataFrames containing the equivalent statistics for each block;
#' and \code{block.ids}, a vector containing the identities of the unique blocks in the same order as \code{per.block}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{model_gene_variances} function in \url{https://libscran.github.io/scran_variances/}.
#'
#' \code{\link{chooseRnaHvgs.se}}, which computes the variances and trend from a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @examples
#' library(Matrix)
#' x <- abs(rsparsematrix(1000, 100, 0.1) * 10)
#' out <- modelGeneVariances(x)
#' out
#'
#' # Throwing in some blocking.
#' block <- sample(letters[1:4], ncol(x), replace=TRUE)
#' out <- modelGeneVariances(x, block=block)
#' out
#'
#' @export
#' @importFrom beachmat initializeCpp
#' @importFrom S4Vectors DataFrame
modelGeneVariances <- function(
    x,
    block=NULL,
    block.average.policy=c("mean", "quantile"),
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    block.quantile=0.5,
    mean.filter=TRUE,
    min.mean=0.1, 
    transform=TRUE, 
    span=0.3,
    use.min.width=FALSE,
    min.width=1,
    min.window.count=200,
    num.threads=1
) {
    block <- .transformFactor(block)

    stats <- model_gene_variances(
        initializeCpp(x, .check.na=FALSE),
        block=block$index,
        nblocks=length(block$names),
        block_average_policy=match.arg(block.average.policy),
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        block_quantile=block.quantile,
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
        statistics = DataFrame(
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
            pb[[i]] <- DataFrame( 
                means = curstats$means,
                variances = curstats$variances,
                fitted = curstats$fitted,
                residuals = curstats$residuals,
                row.names = rownames(x)
            )
        }
        names(pb) <- block$names
        output$per.block <- pb
        output$block.ids <- block$names
    }

    output
}
