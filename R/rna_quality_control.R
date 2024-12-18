#' Quality control for RNA count data
#'
#' Compute per-cell QC metrics from an initialized matrix of RNA counts,
#' and use the metrics to suggest filter thresholds to retain high-quality cells.
#' 
#' @param x A matrix-like object where rows are genes and columns are cells.
#' Values are expected to be counts.
#' @param subsets List of vectors specifying gene subsets of interest, typically for control-like features like mitochondrial genes or spike-in transcripts.
#' Each vector may be logical (whether to keep each row), integer (row indices) or character (row names).
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param metrics List with the same structure as produced by \code{computeRnaQcMetrics}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{metrics}.
#' Alternatively \code{NULL} if all cells are from the same block.
#'
#' For \code{filterRnaQcMetrics}, a blocking factor should be provided if \code{block} was used to construct \code{thresholds}. 
#' @param num.mads Number of median from the median, to define the threshold for outliers in each metric.
#' @param thresholds List with the same structure as produced by \code{suggestRnaQcThresholds}.
#'
#' @return For \code{computeRnaQcMetrics}, a list is returned containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total RNA count for each cell.
#' \item \code{detected}, an integer vector containing the number of detected genes per cell.
#' \item \code{subsets}, a list of numeric vectors containing the proportion of counts in each feature subset.
#' }
#' Each vector is of length equal to the number of cells.
#'
#' For \code{suggestRnaQcThresholds} with \code{block!=NULL}, a list is returned containing:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the lower bound on the sum for each blocking level.
#' \item \code{detected}, a numeric vector containing the lower bound on the number of detected genes for each blocking level.
#' \item \code{subsets}, a list of numeric vectors containing the upper bound on the sum of counts in each feature subset for each blocking level.
#' }
#' Each vector is of length equal to the number of levels in \code{block}.
#'
#' For \code{suggestRnaQcThresholds} with \code{block=NULL}, a list is returned containing:
#' \itemize{
#' \item \code{sum}, a numeric scalar containing the lower bound on the sum.
#' \item \code{detected}, a numeric scalar containing the lower bound on the number of detected genes. 
#' \item \code{subsets}, a numeric vector containing the upper bound on the sum of counts in each feature subset. 
#' }
#'
#' For \code{filterRnaQcMetrics}, a logical vector of length \code{ncol(x)} is returned indicating which cells are of high quality. 
#'
#' @seealso
#' \url{https://libscran.github.io/scran_qc/}, for the rationale of QC filtering on RNA counts.
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
#' str(qc)
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
    ptr <- initializeCpp(x)

    subsets <- as.list(subsets)
    subsets <- lapply(subsets, .toLogical, n=tatami.dim(ptr)[1], names=rownames(x))

    output <- compute_rna_qc_metrics(ptr, subsets, num_threads=num.threads)
    names(output$subsets) <- names(subsets)
    output
}

#' @export
#' @rdname rna_quality_control
suggestRnaQcThresholds <- function(metrics, block=NULL, num.mads=3) {
    block <- .transformFactor(block)
    thresholds <- suggest_rna_qc_thresholds(metrics, block=block$index, num_mads=num.mads)

    names(thresholds$sum) <- block$names
    names(thresholds$detected) <- block$names
    names(thresholds$subsets) <- names(metrics$subsets)
    for (i in seq_along(metrics$subsets)) {
        names(thresholds$subsets[[i]]) <- block$names
    }

    thresholds
}

#' @export
#' @rdname rna_quality_control
filterRnaQcMetrics <- function(thresholds, metrics, block=NULL) {
    block <- .matchBlock(block, names(thresholds$sum))
    filter_rna_qc_metrics(thresholds, metrics, block=block)
}
