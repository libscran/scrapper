#' Fit a mean-variance trend
#'
#' Fit a trend to the per-gene variances with respect to their means, typically from normalized and log-transformed expression values.
#'
#' @param means Numeric vector containing the mean (log-)expression for each gene.
#' @param variances Numeric vector containing the variance in the (log-)expression for each gene.
#' @param mean.filter Logical scalar indicating whether to filter on the means before trend fitting.
#' The assumption is that there is a bulk of low-abundance genes that are uninteresting and should be removed to avoid skewing the windows of the LOWESS smoother.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used.
#' @param min.mean Numeric scalar specifying the minimum mean of genes to use in trend fitting.
#' Genes with lower means do not participate in the LOWESS fit, to ensure that windows are not skewed towards the majority of low-abundance genes.
#' Instead, the fitted values for these genes are defined by extrapolating the left edge of the fitted trend is extrapolated to the origin.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used;
#' this is chosen based on the typical distribution of means of log-expression values across genes.
#' Ignored if \code{mean.filter=FALSE}.
#' @param transform Logical scalar indicating whether a quarter-root transformation should be applied before trend fitting.
#' This transformation is copied from \code{limma::voom} and shrinks all values towards 1, flattening any sharp gradients in the trend for an easier fit.
#' For example, variances are computed from log-expression values typically exhibit a strong \dQuote{hump} in the mean-variance relationship.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used;
#' @param span Numeric scalar specifying the span of the LOWESS smoother, as a proportion of the total number of points.
#' Larger values improve stability at the cost of sensitivity to changes in low-density regions.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used.
#' Ignored if \code{use.min.width=TRUE}.
#' @param use.min.width Logical scalar indicating whether a minimum width constraint should be applied to the LOWESS smoother.
#' This replaces the proportion-based span for defining each window.
#' Instead, the window for each point must be of a minimum width and is extended until it contains a minimum number of points. 
#' Setting this to `TRUE` ensures that sensitivity is maintained in the trend fit at low-density regions for the distribution of means, e.g., at high abundances.
#' It also avoids overfitting from very small windows in high-density intervals. 
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used.
#' @param min.width Minimum width of the window to use when \code{use.min.width=TRUE}.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used;
#' this is chosen based on the typical range of means in single-cell RNA-seq data.
#' @param min.window.count Minimum number of observations in each window.
#' This ensures that each window contains at least a given number of observations for a stable fit.
#' If the minimum width window contains fewer observations, it is extended using the standard LOWESS logic until the minimum number is achieved.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used.
#' Ignored if \code{use.min.width = FALSE}.
#' @param num.threads Number of threads to use.
#' If \code{NULL}, the default value in \code{\link{fitVarianceTrendDefaults}} is used.
#'
#' @return List containing \code{fitted}, a numeric vector containing the fitted values of the trend for each gene;
#' and \code{residuals}, a numeric vector containing the residuals from the trend.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{modelGeneVariances}}, to compute the means and variances on which the trend is fitted.
#'
#' The \code{fit_variance_trend} function in \url{https://libscran.github.io/scran_variances/}.
#'
#' @examples
#' # Setting up some single-cell-like data.
#' mu <- 2^runif(1000, -10, 10)
#' counts <- matrix(rpois(20 * length(mu), lambda=mu), ncol=20)
#'
#' sf <- centerSizeFactors(colSums(counts))
#' normalized <- normalizeCounts(counts, size.factors=sf)
#' stats <- modelGeneVariances(normalized)
#'
#' out <- fitVarianceTrend(stats$statistics$means, stats$statistics$variances)
#' plot(stats$statistics$means, stats$statistics$variances)
#' curve(approxfun(stats$statistics$means, out$fitted)(x), col="red", add=TRUE)
#'
#' @export
fitVarianceTrend <- function(
    means,
    variances,
    mean.filter = NULL,
    min.mean = NULL,
    transform = NULL,
    span = NULL,
    use.min.width = NULL,
    min.width = NULL,
    min.window.count = NULL,
    num.threads = NULL
) {
    fit_variance_trend(
        means,
        variances,
        mean_filter=mean.filter,
        min_mean=min.mean,
        transform=transform,
        span=span,
        use_min_width=use.min.width,
        min_width=min.width,
        min_window_count=min.window.count,
        num_threads=num.threads
    )
}

#' Default parameters for \code{\link{fitVarianceTrend}}
#' @return Named list containing default values for the various function arguments.
#' @author Aaron Lun
#' @examples
#' fitVarianceTrendDefaults()
#' @export
fitVarianceTrendDefaults <- function() fit_variance_trend_defaults()
