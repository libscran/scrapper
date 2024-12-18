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
#' \item \code{detected}, an integer vector containing the number of detected tags per cell.
#' \item \code{subsets}, a list of numeric vectors containing the total count of each control subset. 
#' }
#' Each vector is of length equal to the number of cells.
#'
#' For \code{suggestAdtQcThresholds} with \code{block!=NULL}, a list is returned containing:
#' \itemize{
#' \item \code{detected}, a numeric vector containing the lower bound on the number of detected tags for each blocking level.
#' \item \code{subsets}, a list of numeric vectors containing the upper bound on the sum of counts in each control subset for each blocking level.
#' }
#' Each vector is of length equal to the number of levels in \code{block}.
#'
#' For \code{suggestAdtQcThresholds} with \code{block=NULL}, a list is returned containing:
#' \itemize{
#' \item \code{detected}, a numeric scalar containing the lower bound on the number of detected tags. 
#' \item \code{subsets}, a numeric vector containing the upper bound on the sum of counts in each control subset. 
#' }
#'
#' For \code{filterAdtQcMetrics}, a logical vector of length \code{ncol(x)} is returned indicating which cells are of high quality. 
#'
#' @seealso
#' \url{https://libscran.github.io/scran_qc/}, for the rationale of QC filtering on ADT counts.
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
    ptr <- initializeCpp(x)

    subsets <- as.list(subsets)
    subsets <- lapply(subsets, .toLogical, n=tatami.dim(ptr)[1], names=rownames(x))

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
    block <- .matchBlock(block, names(thresholds$detected))
    filter_adt_qc_metrics(thresholds, metrics, block=block)
}
