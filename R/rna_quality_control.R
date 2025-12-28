#' Quality control for RNA count data
#'
#' Compute per-cell QC metrics from an initialized matrix of RNA counts,
#' and use the metrics to suggest filter thresholds to retain high-quality cells.
#' 
#' @param x A matrix-like object where rows are genes and columns are cells.
#' Values are expected to be counts.
#' @param subsets Named list of vectors specifying gene subsets of interest, typically for control-like features like mitochondrial genes or spike-in transcripts.
#' Each vector may be logical (whether to keep each row), integer (row indices) or character (row names).
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param metrics \link[S4Vectors]{DataFrame} of per-cell QC metrics.
#' This should have the same structure as the return value of \code{computeRnaQcMetrics}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{metrics}.
#' Alternatively \code{NULL} if all cells are from the same block.
#'
#' For \code{filterRnaQcMetrics}, a blocking factor should be provided if \code{block} was used to construct \code{thresholds}. 
#' @param num.mads Number of median from the median, to define the threshold for outliers in each metric.
#' @param thresholds List with the same structure as produced by \code{suggestRnaQcThresholds}.
#'
#' @return For \code{computeRnaQcMetrics}, a \link[S4Vectors]{DataFrame} is returned with one row per cell in \code{x}.
#' This contains the following columns:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total RNA count for each cell.
#' This represents the efficiency of library preparation and sequencing.
#' Low totals indicate that the library was not successfully captured.
#' \item \code{detected}, an integer vector containing the number of detected genes per cell.
#' This also quantifies library preparation efficiency but with greater focus on capturing transcriptional complexity.
#' \item \code{subsets}, a nested DataFrame where each column corresponds to a feature subset and is a numeric vector containing the proportion of counts in that subset.
#' The exact interpretation of which depends on the nature of the subset.
#' For example, if one subset contains all genes on the mitochondrial chromosome, higher proportions are representative of cell damage;
#' the assumption is that cytoplasmic transcripts leak through tears in the cell membrane while the mitochondria are still trapped inside.
#' The proportion of spike-in transcripts can be interpreted in a similar manner, where the loss of endogenous transcripts results in higher spike-in proportions.
#' }
#' Each vector is of length equal to the number of cells.
#'
#' For \code{suggestRnaQcThresholds}, a named list is returned.
#' \itemize{
#' \item If \code{block=NULL}, the list contains:
#' \itemize{
#' \item \code{sum}, a numeric scalar containing the lower bound on the sum.
#' This is defined as \code{num.mads} MADs below the median of the log-transformed metrics across all cells.
#' \item \code{detected}, a numeric scalar containing the lower bound on the number of detected genes. 
#' This is defined as \code{num.mads} MADs below the median of the log-transformed metrics across all cells.
#' \item \code{subsets}, a numeric vector containing the upper bound on the sum of counts in each feature subset. 
#' This is defined as \code{num.mads} MADs above the median across all cells.
#' }
#' \item Otherwise, if \code{block} is supplied, the list contains:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the lower bound on the sum for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' \item \code{detected}, a numeric vector containing the lower bound on the number of detected genes for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' \item \code{subsets}, a list of numeric vectors containing the upper bound on the sum of counts in each feature subset for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' \item \code{block.levels}, a vector containing the unique blocking levels.
#' }
#' Each vector is of length equal to the number of levels in \code{block} and is named accordingly.
#' }
#'
#' For \code{filterRnaQcMetrics}, a logical vector of length \code{ncol(x)} is returned indicating which cells are of high quality. 
#' High-quality cells are defined as those with sums and detected genes above their respective thresholds and subset proportions below the \code{subsets} threshold.
#'
#' @seealso
#' The \code{compute_rna_qc_metrics}, \code{compute_rna_qc_filters} and \code{compute_rna_qc_filters_blocked} functions in \url{https://libscran.github.io/scran_qc/}.
#'
#' \code{\link{quickRnaQc.se}}, to run all of the RNA-related QC functions on a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Mocking up a control set.
#' sub <- list(mito=rbinom(nrow(x), 1, 0.1) > 0)
#'
#' qc <- computeRnaQcMetrics(x, sub)
#' qc
#'
#' filt <- suggestRnaQcThresholds(qc)
#' str(filt)
#'
#' keep <- filterRnaQcMetrics(filt, qc)
#' summary(keep)
#'
#' @export
#' @name rna_quality_control
#' @importFrom beachmat initializeCpp tatami.dim
computeRnaQcMetrics <- function(x, subsets, num.threads = 1) {
    stopifnot(length(subsets) == 0 || !is.null(names(subsets)))
    ptr <- initializeCpp(x, .check.na=FALSE)

    subsets <- as.list(subsets)
    subsets <- lapply(subsets, .subsetToLogical, n=tatami.dim(ptr)[1], names=rownames(x))

    output <- compute_rna_qc_metrics(ptr, subsets, num_threads=num.threads)
    .reformatQcMetrics(output, "subsets", subsets, x)
}

#' @export
#' @rdname rna_quality_control
suggestRnaQcThresholds <- function(metrics, block=NULL, num.mads=3) {
    block <- .transformFactor(block)
    metrics <- .simplifyQcMetrics(metrics)
    thresholds <- suggest_rna_qc_thresholds(metrics, block=block$index, num_mads=num.mads)

    names(thresholds$subsets) <- names(metrics$subsets)

    if (!is.null(block$names)) {
        names(thresholds$sum) <- block$names
        names(thresholds$detected) <- block$names
        for (i in seq_along(metrics$subsets)) {
            names(thresholds$subsets[[i]]) <- block$names
        }
        thresholds$block.levels <- block$names # store it separately to preserve any non-character type.
    }

    thresholds
}

#' @export
#' @rdname rna_quality_control
filterRnaQcMetrics <- function(thresholds, metrics, block=NULL) {
    block <- .matchBlockThresholds(block, thresholds$block.levels)
    metrics <- .simplifyQcMetrics(metrics)
    filter_rna_qc_metrics(thresholds, metrics, block=block)
}

#' @importFrom S4Vectors DataFrame I
.reformatQcMetrics <- function(res, subfield, subsets, mat) {
    if (length(res$subsets) > 0) {
        names(res$subsets) <- names(subsets)
        subdf <- DataFrame(res$subsets, row.names=colnames(mat))
    } else {
        subdf <- make_zero_col_DFrame(ncol(mat))
        rownames(subdf) <- colnames(mat)
    }
    res$subsets <- I(subdf)
    DataFrame(res, row.names=colnames(mat))
}

.simplifyQcMetrics <- function(df) {
    metrics <- as.list(df)
    metrics$subsets <- as.list(metrics$subsets)
    metrics
}

.subsetToLogical <- function(x, n, names) {
    if (is.logical(x)) {
        if (length(x) != n) {
            stop("length of a 'subsets' entry must be equal to the number of rows")
        }
    } else if (is.numeric(x)) {
        if (length(x)) {
            if (anyNA(x) || min(x) < 1 || max(x) > n) {
                stop("'subsets' entry contains out-of-range indices")
            }
        }
        tmp <- logical(n)
        tmp[x] <- TRUE
        x <- tmp
    } else if (is.character(x)) {
        if (is.null(names)) {
            stop("no row names available for matching a 'subsets' entry")
        }
        x <- names %in% x
    }
    x
}

.matchBlockThresholds <- function(block, levels) {
    if (is.null(block) && is.null(levels)) {
        return(NULL)
    } else if (is.null(block)) {
        stop("expected 'block=' to be supplied for blocked 'thresholds'")
    } else if (is.null(levels)) {
        stop("'block=' should be set to NULL for unblocked 'thresholds'")
    }

    m <- match(block, levels)
    if (anyNA(m)) {
        stop("entries of 'block' are not present in 'thresholds'")
    }

    m - 1L
}
