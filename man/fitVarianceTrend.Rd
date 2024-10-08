% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitVarianceTrend.R
\name{fitVarianceTrend}
\alias{fitVarianceTrend}
\title{Fit a mean-variance trend}
\usage{
fitVarianceTrend(
  means,
  variances,
  mean.filter = TRUE,
  min.mean = 0.1,
  transform = TRUE,
  span = 0.3,
  use.min.width = FALSE,
  min.width = 1,
  min.window.count = 200,
  num.threads = 1
)
}
\arguments{
\item{means}{Numeric vector containing the mean (log-)expression for each gene.}

\item{variances}{Numeric vector containing the variance in the (log-)expression for each gene.}

\item{mean.filter}{Logical scalar indicating whether to filter on the means before trend fitting.}

\item{min.mean}{Numeric scalar specifying the minimum mean of genes to use in trend fitting.
Only used if \code{mean.filter=TRUE}.}

\item{transform}{Logical scalar indicating whether a quarter-root transformation should be applied before trend fitting.}

\item{span}{Numeric scalar specifying the span of the LOWESS smoother.
Ignored if \code{use.min.width=TRUE}.}

\item{use.min.width}{Logical scalar indicating whether a minimum width constraint should be applied to the LOWESS smoother.
Useful to avoid overfitting in high-density intervals.}

\item{min.width}{Minimum width of the window to use when \code{use.min.width=TRUE}.}

\item{min.window.count}{Minimum number of observations in each window.
Only used if \code{use.min.width=TRUE}.}

\item{num.threads}{Number of threads to use.}
}
\value{
List containing \code{fitted}, the fitted values of the trend for each gene;
and \code{residuals}, the residuals from the trend.
}
\description{
Fit a trend to the per-cell variances with respect to the mean.
}
\examples{
x <- runif(1000)
y <- 2^rnorm(1000)
out <- fitVarianceTrend(x, y)

plot(x, y)
o <- order(x)
lines(x[o], out$fitted[o], col="red")

}
\seealso{
\url{https://libscran.github.io/scran_variances/}, for the underlying implementation.
}
\author{
Aaron Lun
}
