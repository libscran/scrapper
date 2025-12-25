#' Quick quality control for ADT data in a SummarizedExperiment
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from ADT data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' This calls \code{\link{computeAdtQcMetrics}} on an assay in a SummarizedExperiment,
#' followed by \code{\link{suggestAdtQcThresholds}} and \code{\link{filterAdtQcMetrics}} to identify high-quality cells.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to antibody-derived tags (ADTs) and columns correspond to cells.
#' @param subsets List of subsets of control tags, see \code{?\link{computeAdtQcMetrics}} for more details.
#' @param num.threads Number of threads, to pass to \code{\link{computeAdtQcMetrics}}.
#' @param block Block assignment for each cell, to pass to \code{\link{suggestAdtQcThresholds}} and \code{\link{filterAdtQcMetrics}}.
#' @param more.suggest.args Named list of additional arguments to pass to \code{\link{suggestAdtQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the ADT count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing additional outputs like the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res List returned by \code{\link{computeAdtQcMetrics}}.
#' 
#' @return
#' For \code{quickAdtQc.se}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link{computeAdtQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' For \code{formatComputeAdtQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned with the per-cell QC metrics.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- altExp(getTestAdtData.se(), "ADT")
#' sce <- quickAdtQc.se(sce, subsets=list(igg=grepl("IgG", rownames(sce))))
#' colData(sce)[,c("sum", "detected", "igg.sum")]
#' metadata(sce)$qc$thresholds
#' summary(sce$keep)
#' 
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
quickAdtQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    block = NULL,
    more.suggest.args = list(),
    assay.type = "counts",
    output.prefix = NULL, 
    meta.name = "qc",
    flatten = TRUE
) {
    metrics <- computeAdtQcMetrics(SummarizedExperiment::assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)

    thresholds <- .call(
        suggestAdtQcThresholds,
        list(metrics=metrics),
        list(block=block),
        more.suggest.args
    )

    keep <- filterAdtQcMetrics(thresholds, metrics, block=block)

    df <- formatComputeAdtQcMetricsResult(metrics, flatten=flatten)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::colData(x) <- cbind(SummarizedExperiment::colData(x), df)

    if (!is.null(meta.name)) {
        S4Vectors::metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    x
}

#' @export
#' @rdname quickAdtQc.se
formatComputeAdtQcMetricsResult <- function(compute.res, flatten = TRUE) {
    df <- S4Vectors::DataFrame(sum=compute.res$sum, detected=compute.res$detected)

    if (flatten) {
        for (sub in names(compute.res$subsets)) {
            df[[paste0(sub, ".sum")]] <- compute.res$subsets[[sub]]
        }
    } else {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(df))
        for (sub in names(compute.res$subsets)) {
            tmp[[sub]] <- compute.res$subsets[[sub]]
        }
        df[["sum"]] <- tmp
    }

    df
}
