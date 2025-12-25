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
#' @param more.suggest.args Named list of additional arguments to pass to \code{\link{suggestCrisprQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the CRISPR count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing additional outputs like the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported. 
#' @param compute.res List returned by \code{\link[scrapper]{computeCrisprQcMetrics}}.
#' 
#' @return
#' For \code{quickCrisprQc.se}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeCrisprQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
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
quickCrisprQc.se <- function( 
    x,
    num.threads = 1,
    block = NULL,
    more.suggest.args = list(),
    assay.type = "counts",
    output.prefix = NULL,
    meta.name = "qc"
) {
    metrics <- computeCrisprQcMetrics(SummarizedExperiment::assay(x, assay.type, withDimnames=FALSE), num.threads=num.threads)

    thresholds <- .call(
        suggestCrisprQcThresholds,
        list(metrics),
        list(block=block),
        more.suggest.args
    )

    keep <- filterCrisprQcMetrics(thresholds, metrics, block=block)

    df <- formatComputeCrisprQcMetricsResult(metrics)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::colData(x) <- S4Vectors::cbind(SummarizedExperiment::colData(x), df)

    if (!is.null(meta.name)) {
        S4Vectors::metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    x
}

#' @export
#' @rdname quickCrisprQc.se
formatComputeCrisprQcMetricsResult <- function(compute.res) {
    # Maybe we'll add something more later.
    S4Vectors::DataFrame(compute.res)
}
