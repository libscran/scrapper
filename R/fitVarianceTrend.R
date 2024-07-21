#' Fit a mean-variance trend
#'
#' Fit a trend to the per-cell variances with respect to the mean.
#'
#' @param means Numeric vector containing the mean (log-)expression for each gene.
#' @param variances Numeric vector containing the variance in the (log-)expression for each gene.
#' @param mean.filter Logical scalar indicating whether to filter on the means before trend fitting.
#' @param min.mean Numeric scalar specifying the minimum mean of genes to use in trend fitting.
#' Only used if \code{mean.filter=TRUE}.
#' @param transform Logical scalar indicating whether a quarter-root transformation should be applied before trend fitting.
#' @param span Numeric scalar specifying the span of the LOWESS smoother.
#' Ignored if \code{use.min.width=TRUE}.
#' @param use.min.width Logical scalar indicating whether a minimum width constraint should be applied to the LOWESS smoother.
#' Useful to avoid overfitting in high-density intervals.
#' @param min.width Minimum width of the window to use when \code{use.min.width=TRUE}.
#' @param min.window.count Minimum number of observations in each window.
#' Only used if \code{use.min.width=TRUE}.
#' @param num.threads Number of threads to use.
#'
#' @return List containing \code{fitted}, the fitted values of the trend for each gene;
#' and \code{residuals}, the residuals from the trend.
#'
#' @author Aaron Lun
#'
#' @seealso
#' The \code{fit_variance_trend} function in \url{https://libscran.github.io/scran_variances}.
#'
#' @examples
#' x <- runif(1000)
#' y <- 2^rnorm(1000)
#' out <- fitVarianceTrend(x, y)
#'
#' plot(x, y)
#' o <- order(x)
#' lines(x[o], out$fitted[o], col="red")
#'
#' @export
fitVarianceTrend <- function(means, variances, mean.filter=TRUE, min.mean=0.1, transform=TRUE, span=0.3, use.min.width=FALSE, min.width=1, min.window.count=200, num.threads=1) {
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
