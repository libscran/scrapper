#' Quality control for ADT count data
#'
#' Compute per-cell QC metrics from an initialized matrix of ADT counts,
#' and use the metrics to suggest filter thresholds to retain high-quality cells.
#' 
#' @param x A matrix-like object where rows are ADTs and columns are cells.
#' Values are expected to be counts.
#' @param subsets List of vectors specifying tag subsets of interest, typically control tags like IgGs.
#' Each vector may be logical (whether to keep each row), integer (row indices) or character (row names).
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param metrics List with the same structure as produced by \code{computeAdtQcMetrics}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{metrics}.
#' Alternatively \code{NULL} if all cells are from the same block.
#'
#' For \code{filterAdtQcMetrics}, a blocking factor should be provided if \code{block} was used to construct \code{thresholds}. 
#' @param min.detected.drop Minimum drop in the number of detected features from the median, in order to consider a cell to be of low quality.
#' @param num.mads Number of median from the median, to define the threshold for outliers in each metric.
#' @param thresholds List with the same structure as produced by \code{suggestAdtQcThresholds}.
#'
#' @return For \code{computeAdtQcMetrics}, a list is returned containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total ADT count for each cell.
#' In theory, this represents the efficiency of library preparation and sequencing.
#' Compared to RNA, the sum is less useful as a QC metric for ADT data as it is strongly influenced by biological variation in the abundance of the targeted features.
#' Nonetheless, we compute it for diagnostic purposes.
#' \item \code{detected}, an integer vector containing the number of detected tags per cell.
#' Even though ADTs are typically used in situations where few features are highly abundant (e.g., cell type-specific markers), 
#' we still expect detectable coverage of most features due to ambient contamination, non-specific binding or some background expression.
#' Low numbers of detected tags indicates that library preparation or sequencing depth was suboptimal.
#' \item \code{subsets}, a list of numeric vectors containing the total count of each control subset. 
#' The exact interpretation depends on the nature of the feature subset but the most common use case involves isotype control (IgG) features.
#' IgG antibodies should not bind to anything so a high subset sum suggests that non-specific binding is a problem, e.g., due to antibody conjugates.
#' (Unlike RNA quality control, we do not use proportions here as it is entirely possible for a cell to have low counts for other tags due to the absence of their targeted features;
#' this would result in a high proportion even if the cell has a "normal" level of non-specific binding.)
#' }
#' Each vector is of length equal to the number of cells.
#'
#' For \code{suggestAdtQcThresholds}, a named list is returned:
#' \itemize{
#' \item If \code{block=NULL}, the list contains:
#' \itemize{
#' \item \code{detected}, a numeric scalar containing the lower bound on the number of detected tags. 
#' This is defined as the lower of (i) \code{num.mads} MADs below the median for the log-transformed values across all cells,
#' and (ii) the product of \code{1 - min.detected.drop} and the median across all cells.
#' The latter avoids overly aggressive filtering when the MAD is zero.
#' \item \code{subsets}, a numeric vector containing the upper bound on the sum of counts in each control subset. 
#' This is defined as \code{num.mads} MADs above the median of the log-transformed metrics across all cells.
#' }
#' \item Otherwise, if \code{block} is supplied, the list contains:
#' \itemize{
#' \item \code{detected}, a numeric vector containing the lower bound on the number of detected tags for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' \item \code{subsets}, a list of numeric vectors containing the upper bound on the sum of counts in each control subset for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' }
#' Each vector is of length equal to the number of levels in \code{block} and is named accordingly.
#' }
#'
#' For \code{filterAdtQcMetrics}, a logical vector of length \code{ncol(x)} is returned indicating which cells are of high quality. 
#' High-quality cells are defined as those with numbers of detected tags above the \code{detected} threshold and control subset sums below the \code{subsets} threshold.
#' 
#' @seealso
#' The \code{compute_adt_qc_metrics}, \code{compute_adt_qc_filters} and \code{compute_adt_qc_filters_blocked} functions in \url{https://libscran.github.io/scran_qc/}.
#'
#' \code{\link{quickAdtQc.se}}, to run all of the ADT-related QC functions on a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Mocking up a control set.
#' sub <- list(IgG=rbinom(nrow(x), 1, 0.1) > 0)
#'
#' qc <- computeAdtQcMetrics(x, sub)
#' str(qc)
#'
#' filt <- suggestAdtQcThresholds(qc)
#' str(filt)
#'
#' keep <- filterAdtQcMetrics(filt, qc)
#' summary(keep)
#'
#' @export
#' @name adt_quality_control
#' @importFrom beachmat initializeCpp tatami.dim
computeAdtQcMetrics <- function(x, subsets, num.threads = 1) {
    ptr <- initializeCpp(x, .check.na=FALSE)

    subsets <- as.list(subsets)
    subsets <- lapply(subsets, .subsetToLogical, n=tatami.dim(ptr)[1], names=rownames(x))

    output <- compute_adt_qc_metrics(ptr, subsets, num_threads=num.threads)
    names(output$subsets) <- names(subsets)
    output
}

#' @export
#' @rdname adt_quality_control
suggestAdtQcThresholds <- function(metrics, block=NULL, min.detected.drop=0.1, num.mads=3) {
    block <- .transformFactor(block) 
    thresholds <- suggest_adt_qc_thresholds(metrics, block=block$index, min_detected_drop=min.detected.drop, num_mads=num.mads)

    names(thresholds$detected) <- block$names
    names(thresholds$subsets) <- names(metrics$subsets)
    for (i in seq_along(metrics$subsets)) {
        names(thresholds$subsets[[i]]) <- block$names
    }

    thresholds
}

#' @export
#' @rdname adt_quality_control
filterAdtQcMetrics <- function(thresholds, metrics, block=NULL) {
    block <- .matchBlockThresholds(block, names(thresholds$detected))
    filter_adt_qc_metrics(thresholds, metrics, block=block)
}
