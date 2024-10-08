.transformFactor <- function(f) {
    if (is.null(f)) {
        list(index=NULL, names=NULL)
    } else {
        f <- factor(f)
        list(index=as.integer(f) - 1L, names=levels(f))
    }
}

.toLogical <- function(x, n, names) {
    if (is.logical(x)) {
        stopifnot(identical(length(x), n))
    } else if (is.numeric(x)) {
        stopifnot(min(x) >= 1 && max(x) <= n)
        tmp <- logical(n)
        tmp[x] <- TRUE
        x <- tmp
    } else if (is.character(x)) {
        stopifnot(!is.null(names))
        x <- names %in% x
    }
    x
}

.checkIndices <- function(index) {
    stopifnot(length(dim(index)) == 2L)
    r <- range(index)
    if (!is.finite(r[1]) || r[1] < 1L) {
        stop("'index' should contain finite positive integers")
    }
    if (!is.finite(r[2]) || r[2] > ncol(index)) {
        stop("'index' should contain finite integers no greater than 'ncol(index)'")
    }
}
