#' Delayed log-normalization of a matrix
#'
#' Delayed calculation of log-normalized expression values, typically returned by \code{\link{normalizeCounts}}.
#'
#' @param x Count matrix to be normalized. 
#' @param size.factors Numeric vector of size factors, of length equal to the number of columns of \code{x}.
#' @param pseudo.count Number specifying the pseudo-count to add prior to log-transformation.
#' @param log.base Number specifying the base of the log-transformation. 
#' 
#' @return An instance of a LogNormalizedMatrix(Seed).
#'
#' @details
#' This is based on the \pkg{DelayedArray} framework and 
#'
#' @author Aaron Lun
#' @examples
#' mat <- matrix(rpois(1000, lambda=2), ncol=10)
#' sf <- centerSizeFactors(colSums(mat))
#' norm <- LogNormalizedMatrix(mat, sf)
#' norm
#'
#' # Also works with sparse matrices.
#' library(Matrix)
#' smat <- abs(rsparsematrix(50, 20, density=0.21))
#' ssf <- centerSizeFactors(colSums(smat))
#' snorm <- LogNormalizedMatrix(smat, ssf)
#' snorm
#'
#' @aliases
#' LogNormalizedMatrix-class
#' LogNormalizedMatrixSeed-class
#' dim,LogNormalizedMatrixSeed-method
#' dimnames,LogNormalizedMatrixSeed-method
#' type,LogNormalizedMatrixSeed-method
#' is_sparse,LogNormalizedMatrixSeed-method
#' extract_array,LogNormalizedMatrixSeed-method
#' extract_sparse_array,LogNormalizedMatrixSeed-method
#' matrixClass,LogNormalizedMatrixSeed-method
#' DelayedArray,LogNormalizedMatrixSeed-method
#' initializeCpp,LogNormalizedMatrixSeed-method
#'
#' @export
LogNormalizedMatrix <- function(x, size.factors, pseudo.count=1, log.base=2) {
    DelayedArray(LogNormalizedMatrixSeed(x, size.factors=size.factors, pseudo.count=pseudo.count, log.base=log.base))
}

#' @export
#' @rdname LogNormalizedMatrix
#' @importFrom methods new
LogNormalizedMatrixSeed <- function(x, size.factors, pseudo.count=1, log.base=2) {
    new("LogNormalizedMatrixSeed",
        seed=x,
        size.factors=as.double(size.factors), 
        pseudo.count=as.double(pseudo.count),
        log.base=as.double(log.base)
    )
}

#' @export
#' @import methods
setClass("LogNormalizedMatrixSeed", slots=c(seed="ANY", size.factors="numeric", pseudo.count="numeric", log.base="numeric"))

#' @importFrom S4Vectors setValidity2
setValidity2("LogNormalizedMatrixSeed", function(object) {
    if (ncol(object@seed) != length(object@size.factors)) {
        return("length of 'size.factors' should equal number of columns")
    }
})

#' @export
setMethod("dim", "LogNormalizedMatrixSeed", function(x) dim(x@seed))

#' @export
setMethod("dimnames", "LogNormalizedMatrixSeed", function(x) dimnames(x@seed))

#' @export
#' @importFrom DelayedArray type
setMethod("type", "LogNormalizedMatrixSeed", function(x) "double")

#' @export
#' @importFrom SparseArray is_sparse
setMethod("is_sparse", "LogNormalizedMatrixSeed", function(x) is_sparse(x@seed) && x@pseudo.count == 1)

.slice_size_factors <- function(x, index) { 
    sf <- x@size.factors
    if (!is.null(index[[2]])) {
        sf <- sf[index[[2]]]
    }
    sf
}

#' @importFrom Matrix t
.log_normalize <- function(mat, sf, pseudo.count, log.base) {
    norm <- t(t(mat)/sf)
    if (pseudo.count == 1) {
        log1p(norm) / log(log.base)
    } else {
        log(norm + pseudo.count, log.base)
    }
}

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "LogNormalizedMatrixSeed", function(x, index) {
    core <- extract_array(x@seed, index)
    .log_normalize(core, .slice_size_factors(x, index), pseudo.count=x@pseudo.count, log.base=x@log.base)
})

#' @export
#' @importFrom SparseArray extract_sparse_array
setMethod("extract_sparse_array", "LogNormalizedMatrixSeed", function(x, index) {
    core <- extract_sparse_array(x@seed, index)
    .log_normalize(core, .slice_size_factors(x, index), pseudo.count=x@pseudo.count, log.base=x@log.base)
})

#' @export
setClass("LogNormalizedMatrix", contains="DelayedMatrix", slots=c(seed="LogNormalizedMatrixSeed"))

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "LogNormalizedMatrixSeed", function(seed) new_DelayedArray(seed, Class = "LogNormalizedMatrix"))

#' @export
#' @importFrom beachmat initializeCpp
setMethod("initializeCpp", "LogNormalizedMatrixSeed", function(x, ...) {
    seed <- initializeCpp(x@seed)
    initialize_LogNormalizedMatrix(seed, x@size.factors, x@pseudo.count, x@log.base)
})
