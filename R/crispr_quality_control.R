#' Quality control for CRISPR count data
#'
#' Compute per-cell QC metrics from an initialized matrix of CRISPR counts,
#' and use the metrics to suggest filter thresholds to retain high-quality cells.
#' 
#' @param x A matrix-like object where rows are CRISPRs and columns are cells.
#' Values are expected to be counts.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param metrics \link[S4Vectors]{DataFrame} of per-cell QC metrics.
#' This should have the same structure as the return value of \code{computeCrisprQcMetrics}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{metrics}.
#' Alternatively \code{NULL} if all cells are from the same block.
#'
#' For \code{filterCrisprQcMetrics}, a blocking factor should be provided if \code{block} was used to construct \code{thresholds}. 
#' @param num.mads Number of median from the median, to define the threshold for outliers in each metric.
#' @param thresholds List with the same structure as produced by \code{suggestCrisprQcThresholds}.
#'
#' @return For \code{computeCrisprQcMetrics}, a \link[S4Vectors]{DataFrame} is returned with one row per cell in \code{x}.
#' This contains the following columns:
#' \itemize{
#' \item \code{sum}, a numeric vector containing the total CRISPR count for each cell.
#' Low counts indicate that the cell was not successfully transfected with a construct or that library preparation and sequencing failed.
#' \item \code{detected}, an integer vector containing the number of detected guides per cell.
#' In theory, this should be 1, as each cell should express no more than one guide construct.
#' However, ambient contamination may introduce non-zero counts for multiple guides, without necessarily interfering with downstream analyses.
#' As such, this metric is less useful for guide data, though we compute it anyway.
#' \item \code{max.value}, a numeric vector containing the count for the most abundant guide in cell.
#' Low values indicate that the cell was not successfully transfected or that library preparation and sequencing failed.
#' \item \code{max.index}, an integer vector containing the row index of the most abundant guide in cell.
#' }
#' Each vector is of length equal to the number of cells.
#'
#' For \code{suggestCrisprQcThresholds}, a named list is returned.
#' \itemize{
#' \item If \code{block=NULL}, the list contains:
#' \itemize{
#' \item \code{max.value}, a numeric scalar containing the lower bound on the maximum count. 
#' This is defined as \code{num.mads} MADs below the median of the log-transformed metrics across cells with high maximum proportions (see Details).
#' }
#' \item Otherwise, if \code{block} is supplied, the list contains:
#' \itemize{
#' \item \code{max.value}, a numeric vector containing the lower bound on the maximum counts for each blocking level.
#' Here, the threshold is computed independently for each block, using the same method as the unblocked case.
#' \item \code{block.ids}, a vector containing the identities of the unique blocks.
#' }
#' Each vector is of length equal to the number of levels in \code{block} and is named accordingly.
#' }
#'
#' For \code{filterCrisprQcMetrics}, a logical vector of length \code{ncol(x)} is returned indicating which cells are of high quality. 
#' High-quality cells are defined as those with maximum counts above the \code{max.value} threshold.
#'
#' @details
#' In CRISPR data, a cell is considered to be of low quality if it has a low count for its most abundant guide.
#' However, directly defining a MAD-based outlier threshold on the maximum count is somewhat tricky as unsuccessful transfection can be common.
#' This often results in a large subpopulation with low maximum counts, inflating the MAD and compromising the threshold calculation.
#' Instead, we use the following approach:
#' \itemize{
#' \item Compute the proportion of counts in the most abundant guide (i.e., the maximum proportion) in each cell.
#' Cells that were successfully transfected should have high maximum proportions.
#' In contrast, unsuccessfully transfected cells will be dominated by ambient contamination and have low proportions.
#' \item Subset the dataset to only retain those cells with maximum proportions above the median.
#' This assumes that at least 50% of cells were successfully transfected.
#' Thus, we remove all of the unsucessful transfections and enrich for mostly-high-quality cells.
#' \item Define a MAD-based threshold for low outliers on the log-transformed maximum count within the subset (see `choose_filter_thresholds()` for details).
#' This is now possible as we can assume that most of the remaining cells are of high quality.
#' }
#' Note that the maximum proportion is only used to define the subset for threshold calculation.
#' Once the maximum count threshold is computed, it is applied to all cells regardless of their maximum proportions.
#' This ensures that we correctly remove cells with low coverage, even if the proportion is high.
#' It also allows us to retain cells transfected with multiple guides, as long as the maximum is high enough -
#' such cells are not necessarily uninteresting, e.g., for examining interaction effects, so we will err on the side of caution and leave them in.
#'
#' @seealso
#' The \code{compute_crispr_qc_metrics}, \code{compute_crispr_qc_filters} and \code{compute_crispr_qc_filters_blocked} functions in \url{https://libscran.github.io/scran_qc/}.
#'
#' \code{\link{quickCrisprQc.se}}, to run all of the CRISPR-related QC functions on a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(100, 100, 0.1) * 100))
#'
#' qc <- computeCrisprQcMetrics(x)
#' qc
#'
#' filt <- suggestCrisprQcThresholds(qc)
#' str(filt)
#'
#' keep <- filterCrisprQcMetrics(filt, qc)
#' summary(keep)
#'
#' @export
#' @name crispr_quality_control
#' @importFrom S4Vectors DataFrame
computeCrisprQcMetrics <- function(x, num.threads = 1) {
    y <- initializeCpp(x, .check.na=FALSE)
    output <- compute_crispr_qc_metrics(y, num_threads=num.threads)
    output$max.index <- output$max.index + 1L
    DataFrame(output, row.names=colnames(x))
}

#' @export
#' @rdname crispr_quality_control
suggestCrisprQcThresholds <- function(metrics, block=NULL, num.mads=3) {
    block <- .transformFactor(block)
    metrics <- as.list(metrics)
    metrics$max.index <- metrics$max.index - 1L # restore 0-based indexing.
    thresholds <- suggest_crispr_qc_thresholds(metrics, block=block$index, num_mads=num.mads)

    if (!is.null(block$names)) {
        names(thresholds$max.value) <- block$names
        thresholds$block.ids <- block$names # store it separately to preserve any non-character type.
    }

    thresholds
}

#' @export
#' @rdname crispr_quality_control
filterCrisprQcMetrics <- function(thresholds, metrics, block=NULL) {
    block <- .matchBlockThresholds(block, thresholds$block.ids)
    metrics <- as.list(metrics)
    metrics$max.index <- metrics$max.index - 1L # restore 0-based indexing.
    .checkThresholdNames(thresholds, c("max.value"), NULL, NULL)
    filter_crispr_qc_metrics(thresholds, metrics, block=block)
}
