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
#' @param thresholds List containing pre-defined thresholds for each QC metric.
#' This should contain zero, one or more of:
#' \itemize{
#' \item \code{detected}, a number containing the lower bound on the number of detected tags. 
#' \item \code{subset.sum}, a named list of numbers containing the upper bound on the sum of counts in each control subset. 
#' }
#' If any element is missing, it will be filled in using the corresponding value from \code{\link{suggestAdtQcThresholds}}.
#'
#' If \code{block} is provided, each numeric scalar may instead be a named numeric vector containing the threshold to be used in each block.
#' Otherwise, the threshold in the numeric scalar will be applied to all blocks.
#'
#' If \code{NULL}, all thresholds are computed by \code{suggestAdtQcThresholds}.
#'
#' If \code{thresholds} is wrapped in an \link{AsIs}, it will be used verbatim, and \code{suggestAdtQcThresholds} will not be called.
#' @param more.suggest.args Named list of additional arguments to pass to \code{\link{suggestAdtQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the ADT count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{\link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing additional outputs like the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{\link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res \link[S4Vectors]{DataFrame} returned by \code{\link{computeAdtQcMetrics}}.
#' @param filter.cells Boolean indicating whether to filter the output to retain only high-quality cells.
#' 
#' @return
#' For \code{quickAdtQc.se}, a copy of \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link{computeAdtQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' If \code{filter.cells=TRUE}, the output object is filtered so that it only contains columns for which \code{keep} is true.
#' By default, no filtering is performed so the number and order of the columns in the output is the same as \code{x}.
#'
#' For \code{formatComputeAdtQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned with the per-cell QC metrics.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- altExp(getTestAdtData.se(), "ADT")
#' sce <- quickAdtQc.se(sce, subsets=list(igg=grepl("IgG", rownames(sce))))
#' colData(sce)[,c("sum", "detected", "subset.sum.igg")]
#' metadata(sce)$qc$thresholds
#' summary(sce$keep)
#' 
#' @export
#' @importFrom S4Vectors cbind metadata metadata<-
quickAdtQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    thresholds = NULL,
    block = NULL,
    more.suggest.args = list(),
    assay.type = "counts",
    output.prefix = NULL, 
    meta.name = "qc",
    flatten = TRUE,
    filter.cells = FALSE
) {
    metrics <- computeAdtQcMetrics(SummarizedExperiment::assay(x, assay.type), subsets, num.threads=num.threads)

    if (is(thresholds, "AsIs")) {
        thresholds <- unclass(thresholds)
    } else {
        thresholds.0 <- .call(
            suggestAdtQcThresholds,
            list(metrics=metrics),
            list(block=block),
            more.suggest.args
        )
        if (is.null(thresholds)) {
            thresholds <- thresholds.0
        } else {
            names(thresholds)[names(thresholds) == "subset.sum"] <- "subsets"
            thresholds <- .substituteThresholds(thresholds.0, thresholds)
        }
    }

    keep <- filterAdtQcMetrics(thresholds, metrics, block=block)

    df <- formatComputeAdtQcMetricsResult(metrics, flatten=flatten)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::colData(x) <- .cbindOverwrite(SummarizedExperiment::colData(x), df)

    if (!is.null(meta.name)) {
        names(thresholds)[names(thresholds) == "subsets"] <- "subset.sum"
        metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    if (filter.cells) {
        x <- .delayifyAssays(x)
        x <- x[,keep]
    }

    x
}

#' @export
#' @rdname quickAdtQc.se
#' @importFrom S4Vectors DataFrame cbind
formatComputeAdtQcMetricsResult <- function(compute.res, flatten = TRUE) {
    if (flatten) {
        extra <- compute.res$subsets
        colnames(extra) <- sprintf("subset.sum.%s", colnames(extra))
        cbind(DataFrame(sum=compute.res$sum, detected=compute.res$detected), extra)
    } else {
        colnames(compute.res)[colnames(compute.res) == "subsets"] <- "subset.sum"
        compute.res
    }
}
