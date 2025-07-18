% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correctMnn.R
\name{correctMnn}
\alias{correctMnn}
\title{Batch correction with mutual nearest neighbors}
\usage{
correctMnn(
  x,
  block,
  num.neighbors = 15,
  num.steps = 1,
  merge.policy = c("rss", "size", "variance", "input"),
  num.mads = NULL,
  robust.iterations = NULL,
  robust.trim = NULL,
  mass.cap = NULL,
  order = NULL,
  reference.policy = NULL,
  BNPARAM = AnnoyParam(),
  num.threads = 1
)
}
\arguments{
\item{x}{Numeric matrix where rows are dimensions and columns are cells,
typically containing low-dimensional coordinates (e.g., from \code{\link{runPca}}).}

\item{block}{Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.}

\item{num.neighbors}{Integer scalar specifying the number of neighbors to use when identifying MNN pairs.}

\item{num.steps}{Integer scalar specifying the number of steps for the recursive neighbor search to compute the center of mass.}

\item{merge.policy}{String specifying the policy to use to choose the order of batches to merge.
This can be based on the size of the batch (\code{"size"}),
the variance within each batch (\code{"variance"}), 
the residual sum of squares of each batch (\code{"rss"}),
or the input order (\code{"input"}).}

\item{num.mads}{Deprecated and ignored.}

\item{robust.iterations}{Deprecated and ignored.}

\item{robust.trim}{Deprecated and ignored.}

\item{mass.cap}{Deprecated and ignored.}

\item{order}{Deprecated and ignored, the merge order is now always automatically determined.}

\item{reference.policy}{Deprecated, use \code{merge.policy} instead.}

\item{BNPARAM}{A \link[BiocNeighbors]{BiocNeighborParam} object specifying the nearest-neighbor algorithm to use.}

\item{num.threads}{Integer scalar specifying the number of threads to use.}
}
\value{
List containing \code{corrected}, a numeric matrix of the same dimensions as \code{x}, containing the corrected values.
}
\description{
Apply mutual nearest neighbor (MNN) correction to remove batch effects from a low-dimensional matrix.
}
\examples{
# Mocking up a dataset with multiple batches.
x <- matrix(rnorm(10000), nrow=10)
b <- sample(3, ncol(x), replace=TRUE)
x[,b==2] <- x[,b==2] + 3
x[,b==3] <- x[,b==3] + 5
lapply(split(colMeans(x), b), mean) # indeed the means differ...

corrected <- correctMnn(x, b)
str(corrected)
lapply(split(colMeans(corrected$corrected), b), mean) # now merged.

}
\seealso{
\url{https://libscran.github.io/mnncorrect/}, for more details on the underlying implementation.
}
\author{
Aaron Lun
}
