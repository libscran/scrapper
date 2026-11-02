#' Quick quality control for CRISPR data in a SummarizedExperiment
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from CRISPR data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' This calls \code{\link{computeCrisprQcMetrics}} on an assay in a SummarizedExperiment,
#' followed by \code{\link{suggestCrisprQcThresholds}} and \code{\link{filterCrisprQcMetrics}} to identify high-quality cells.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to CRISPR guides and columns correspond to cells.
#' @param num.threads Number of threads, to pass to \code{\link{computeCrisprQcMetrics}}.
#' @param block Block assignment for each cell, to pass to \code{\link{suggestCrisprQcThresholds}} and \code{\link{filterCrisprQcMetrics}}.
#' @param thresholds List containing pre-defined thresholds for each QC metric.
#' This should contain zero, one or more of:
#' \itemize{
#' \item \code{max.value}, a number containing the lower bound on the maximum count. 
#' }
#' If any element is missing, it will be filled in using the corresponding value from \code{\link{suggestCrisprQcThresholds}}.
#'
#' If \code{block} is provided, each numeric scalar may instead be a named numeric vector containing the threshold to be used in each block.
#' Otherwise, the threshold in the numeric scalar will be applied to all blocks.
#'
#' If \code{NULL}, all thresholds are computed by \code{suggestCrisprQcThresholds}.
#'
#' If \code{thresholds} is wrapped in an \link{AsIs}, it will be used verbatim, and \code{suggestCrisprQcThresholds} will not be called.
#' @param more.suggest.args Named list of additional arguments to pass to \code{\link{suggestCrisprQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the CRISPR count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{\link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing additional outputs like the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported. 
#' @param compute.res \link[S4Vectors]{DataFrame} returned by \code{\link[scrapper]{computeCrisprQcMetrics}}.
#' @param filter.cells Boolean indicating whether to filter the output to retain only high-quality cells.
#' 
#' @return
#' For \code{quickCrisprQc.se}, a copy of \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeCrisprQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' If \code{filter.cells=TRUE}, the output object is filtered so that it only contains columns for which \code{keep} is true.
#' By default, no filtering is performed so the number and order of the columns in the output is the same as \code{x}.
#'
#' For \code{formatComputeCrisprQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned with the per-cell QC metrics.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- altExp(getTestCrisprData.se(), "CRISPR Guide Capture")
#' sce <- quickCrisprQc.se(sce)
#' colData(sce)[,c("sum", "detected", "max.value", "max.index")]
#' metadata(sce)$qc$thresholds
#' summary(sce$keep)
#' 
#' @export
#' @importFrom S4Vectors cbind metadata metadata<-
quickCrisprQc.se <- function( 
    x,
    num.threads = 1,
    thresholds = NULL,
    block = NULL,
    more.suggest.args = list(),
    assay.type = "counts",
    output.prefix = NULL,
    meta.name = "qc",
    filter.cells = FALSE
) {
    metrics <- computeCrisprQcMetrics(SummarizedExperiment::assay(x, assay.type), num.threads=num.threads)

    if (is(thresholds, "AsIs")) {
        thresholds <- unclass(thresholds)
    } else {
        thresholds.0 <- .call(
            suggestCrisprQcThresholds,
            list(metrics),
            list(block=block),
            more.suggest.args
        )
        if (is.null(thresholds)) {
            thresholds <- thresholds.0
        } else {
            thresholds <- .substituteThresholds(thresholds.0, thresholds)
        }
    }

    keep <- filterCrisprQcMetrics(thresholds, metrics, block=block)

    df <- formatComputeCrisprQcMetricsResult(metrics)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::colData(x) <- .cbindOverwrite(SummarizedExperiment::colData(x), df)

    if (!is.null(meta.name)) {
        metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    if (filter.cells) {
        x <- .delayifyAssays(x)
        x <- x[,keep]
    }

    x
}

#' @export
#' @rdname quickCrisprQc.se
formatComputeCrisprQcMetricsResult <- function(compute.res) {
    # Maybe we'll add something more later.
    compute.res
}
