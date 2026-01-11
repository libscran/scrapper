.transformFactor <- function(f) {
    if (is.null(f)) {
        list(index=NULL, names=NULL)
    } else {
        # Don't use 'factor', to preserve the type of 'f' in the 'names'.
        lev <- sort(unique(f))
        list(index=match(f, lev) - 1L, names=lev)
    }
}

#' @importFrom methods is
.checkSEX <- function(x, alt) {
    if (is(x, "SummarizedExperiment")) {
        stop("SummarizedExperiment inputs are not supported, use '", alt, "()' or extract the relevant 'assay()' instead")
    }
}

.checkNeighborResults <- function(index, distance) {
    stopifnot(length(dim(index)) == 2L)
    r <- range(index)
    if (!is.finite(r[1]) || r[1] < 1L) {
        stop("'index' should contain finite positive integers")
    }
    if (!is.finite(r[2]) || r[2] > ncol(index)) {
        stop("'index' should contain finite integers no greater than 'ncol(index)'")
    }
    if (!is.null(distance) && !identical(dim(index), dim(distance))) {
        stop("'index' and 'distance' should have the same dimensions")
    }
}

# We keep arguments that are explicitly named in 'more.args' for back-compatibility,
# in case some of these are promoted into top-level arguments at a later date.
.call <- function(FUN, essential.args, named.args, more.args) {
    all.args <- c(
        essential.args,
        named.args[!(names(named.args) %in% names(more.args))],
        more.args
    )
    do.call(FUN, all.args)
}

.sanitizeAltexpAssays <- function(altexps, all.altexps, default.assay.type) {
    if (!is.null(names(altexps))) {
        altexps[!duplicated(names(altexps))]
    } else {
        altexps <- unique(altexps)
        output <- rep(default.assay.type, length(altexps))
        if (is.numeric(altexps)) {
            names(output) <- all.altexps[altexps]
        } else {
            names(output) <- altexps
        }
        output
    }
}

#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
.addTransposedReddim <- function(x, name, mat, delayed) {
    if (delayed) {
        mat <- DelayedArray(mat)
    }
    SingleCellExperiment::reducedDim(x, name) <- t(mat)
    x
}

#' @importFrom Matrix t
#' @importFrom methods is
#' @importClassesFrom DelayedArray DelayedArray
.getTransposedReddim <- function(x, name) {
    if (is.null(names(name))) {
        mat <- SingleCellExperiment::reducedDim(x, name, withDimnames=FALSE)
    } else {
        mat <- SingleCellExperiment::reducedDim(SingleCellExperiment::altExp(x, names(name), withDimnames=FALSE), name, withDimnames=FALSE)
    }

    mat <- t(mat)

    if (is.matrix(mat)) {
        colnames(mat) <- colnames(x)
        return(mat)
    } else if (is(mat, "DelayedArray")) {
        # Possibly a no-op if .add_transposed_reddim was set with delayed=TRUE.
        if (is.matrix(mat@seed)) {
            out <- mat@seed
            colnames(out) <- colnames(x)
            return(out)
        }
    }

    colnames(mat) <- colnames(x)
    return(as.matrix(mat))
}
