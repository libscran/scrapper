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
#' This can either be \code{"mean"} (to compute a weighted mean), \code{"quantile"} (to compute a quantile) or \code{"none"} (no averaging).
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Ignored if \code{block = NULL}.
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{block} is not \code{NULL} and \code{block.average.policy = "mean"}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{block} is not \code{NULL}, \code{block.average.policy = "mean"} and \code{block.weight.policy = "variable"}.
#' @param block.quantile Number specifying the probability of the quantile of statistics across blocks. 
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used; this is set to 0.5, i.e., the median of per-block statistics.
#' Only used if \code{block} is not \code{NULL} and \code{block.average.policy = "quantile"}.
#' @param fit.trend Boolean indicating whether a mean-variance trend should be fitted.
#' If \code{FALSE}, only the means and variances are computed.
#' This can occasionally be useful when the trend is computed separately (e.g., spike-ins).
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' @param mean.filter Logical scalar indicating whether to filter on the means before trend fitting.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param min.mean Number specifying the minimum mean of genes to use in trend fitting.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param transform Logical scalar indicating whether a quarter-root transformation should be applied before trend fitting.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' @param span Numeric scalar specifying the span of the LOWESS smoother, as a proportion of the total number of points.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param use.min.width Logical scalar indicating whether a minimum width constraint should be applied to the LOWESS smoother.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param min.width Minimum width of the window to use when \code{use.min.width=TRUE}.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param min.window.count Minimum number of observations in each window.
#' See the argument of the same name in \code{\link{fitVarianceTrend}} for more details.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
#' Only used if \code{fit.trend = TRUE}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' If \code{NULL}, the default value in \code{\link{modelGeneVariancesDefaults}} is used.
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
#' If \code{block} is supplied and \code{block.average.policy = "none"}, the \code{statistics} DataFrame will have no columns.
#'
#' If \code{fit.trend = FALSE}, the \code{fitted} and \code{residuals} columns will be omitted from all DataFrames.
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
    block = NULL,
    block.average.policy = NULL,
    block.weight.policy = NULL,
    variable.block.weight = NULL,
    block.quantile = NULL,
    fit.trend = TRUE,
    mean.filter = TRUE,
    min.mean = 0.1, 
    transform = TRUE, 
    span = 0.3,
    use.min.width = FALSE,
    min.width = 1,
    min.window.count = 200,
    num.threads = 1
) {
    .checkSEX(x, "chooseRnaHvgs.se")
    block <- .transformFactor(block)

    computed <- model_gene_variances(
        initializeCpp(x, .check.na=FALSE),
        block=block$index,
        nblocks=length(block$names),
        block_average_policy = block.average.policy,
        block_weight_policy = block.weight.policy,
        variable_block_weight=variable.block.weight,
        block_quantile=block.quantile,
        fit_trend = fit.trend,
        mean_filter=mean.filter,
        min_mean=min.mean,
        transform=transform,
        span=span,
        use_min_width=use.min.width,
        min_width=min.width,
        min_window_count=min.window.count,
        num_threads=num.threads
    )

    if ("means" %in% names(computed)) {
        stats <- DataFrame(means = computed$means, variances = computed$variances)
        if (fit.trend) {
            stats$fitted <- computed$fitted
            stats$residuals <- computed$residuals
        }
    } else {
        # Otherwise, we're blocking and no average was requested.
        stats <- make_zero_col_DFrame(nrow(x))
    }
    rownames(stats) <- rownames(x)
    output <- list(statistics = stats)

    if (!is.null(block$index)) {
        pb <- computed$per.block
        for (i in seq_along(pb)) {
            curcomp <- pb[[i]]
            curstats <- DataFrame(means = curcomp$means, variances = curcomp$variances, row.names = rownames(x))
            if (fit.trend) {
                curstats$fitted <- curcomp$fitted
                curstats$residuals <- curcomp$residuals
            }
            pb[[i]] <- curstats
        }
        names(pb) <- block$names
        output$per.block <- pb
        output$block.ids <- block$names
    }

    output
}

#' Default parameters for \code{\link{modelGeneVariances}}
#' @return Named list of default values for various function arguments. 
#' @author Aaron Lun
#' @examples
#' modelGeneVariancesDefaults()
#' @export
modelGeneVariancesDefaults <- function() model_gene_variances_defaults()
