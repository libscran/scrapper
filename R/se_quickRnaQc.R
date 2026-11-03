#' Quick quality control for RNA data in a SummarizedExperiment
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from RNA data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' This calls \code{\link{computeRnaQcMetrics}} on an assay in a SummarizedExperiment,
#' followed by \code{\link{suggestRnaQcThresholds}} and \code{\link{filterRnaQcMetrics}} to identify high-quality cells.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param subsets List of subsets of control genes, see \code{?\link{computeRnaQcMetrics}} for more details.
#' @param num.threads Number of threads, to pass to \code{\link{computeRnaQcMetrics}}.
#' @param block Block assignment for each cell, to pass to \code{\link{suggestRnaQcThresholds}} and \code{\link{filterRnaQcMetrics}}.
#' @param thresholds List containing pre-defined thresholds for each QC metric.
#' This should contain zero, one or more of:
#' \itemize{
#' \item \code{sum}, a number containing the lower bound on the sum.
#' \item \code{detected}, a number containing the lower bound on the number of detected genes. 
#' \item \code{subset.proportion}, a named list of numbers containing the upper bound on the proportion of counts in each feature subset. 
#' }
#' If any element is missing, it will be filled in using the corresponding value from \code{\link{suggestRnaQcThresholds}}.
#'
#' If \code{block} is provided, each number may instead be a named numeric vector containing the threshold to be used in each block.
#' Otherwise, the threshold in the number will be applied to all blocks.
#'
#' If \code{NULL}, all thresholds are computed by \code{suggestRnaQcThresholds}.
#'
#' If \code{thresholds} is wrapped in an \link{AsIs}, it will be used verbatim, and \code{suggestRnaQcThresholds} will not be called.
#' @param more.suggest.args Named list of additional arguments to pass to \code{\link{suggestRnaQcThresholds}}.
#' @param altexp.proportions Alternative experiments for which to compute QC metrics.
#' This is typically used to refer to alternative experiments holding spike-in data.
#' For each alternative experiment, the proportion is defined as \eqn{X/(X+Y)} where \eqn{X} is the alternative experiment's total and \eqn{Y} is the RNA total.
#' These proportions will be used for filtering in the same manner as the proportions computed from \code{subsets}.
#'
#' More specifically, \code{altexp.proportions} should be an unnamed integer or character vector containing the names/indices of the alternative experiments of interest.
#' The assay to use from each alternative experiment is determined by \code{assay.type}.
#'
#' Alternatively, \code{altexp.proportions} may be a named integer or character vector.
#' Each name specifies an alternative experiment while each value is the index/name of the assay to use from that experiment.
#' 
#' Only relevant if \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the RNA count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{\link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing the additional outputs such as the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported. 
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{\link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res \link[S4Vectors]{DataFrame} returned by \code{\link{computeRnaQcMetrics}}.
#' @param filter.cells Boolean indicating whether to filter the output to retain only high-quality cells.
#' 
#' @return
#' For \code{quickRnaQc.se}, a copy of \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link{computeRnaQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' If \code{filter.cells=TRUE}, the output object is filtered so that it only contains columns for which \code{keep} is true.
#' By default, no filtering is performed so the number and order of the columns in the output is the same as \code{x}.
#'
#' If \code{altexp.proportions} is specified, QC metrics are added to the \code{colData} of the specified alternative experiments in the output object.
#'
#' For \code{computeRnaQcMetricsWithAltExps}, a list is returned containing:
#' \itemize{
#' \item \code{main}, the result of calling \code{\link{computeRnaQcMetrics}} on the RNA count matrix in \code{x}.
#' The proportion of counts in each alternative experiment is added to the \code{subsets}.
#' \item \code{altexp}, a named list of length equal to \code{altexp.proportions}.
#' Each inner list is the result of calling \code{computeRnaQcMetrics} on the RNA count matrix of the corresponding alternative experiment of \code{x}.
#' }
#'
#' For \code{formatComputeRnaQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned containing the per-cell QC metrics.
#'
#' @author Aaron Lun
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se()
#' sce <- quickRnaQc.se(sce, subsets=list(mito=grepl("^mt", rownames(sce))))
#' colData(sce)[,c("sum", "detected", "subset.proportion.mito")]
#' metadata(sce)$qc$thresholds
#' summary(sce$keep)
#'
#' # Computing spike-in proportions, if available.
#' sce <- getTestRnaData.se()
#' sce <- quickRnaQc.se(
#'    sce,
#'    subsets=list(mito=grepl("^mt", rownames(sce))),
#'    altexp.proportions="ERCC"
#' )
#' colData(sce)[,c("sum", "detected", "subset.proportion.mito", "subset.proportion.ERCC")]
#' colData(altExp(sce, "ERCC"))[,c("sum", "detected")]
#'
#' @export
#' @importFrom S4Vectors cbind metadata metadata<-
quickRnaQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    thresholds = NULL,
    block = NULL,
    more.suggest.args = list(),
    altexp.proportions = NULL,
    assay.type = "counts",
    output.prefix = NULL, 
    meta.name = "qc",
    flatten = TRUE,
    filter.cells = FALSE
) {
    metrics <- computeRnaQcMetricsWithAltExps(
        x,
        subsets,
        altexp.proportions=altexp.proportions,
        assay.type=assay.type,
        num.threads=num.threads
    )

    if (is(thresholds, "AsIs")) {
        thresholds <- unclass(thresholds)
    } else {
        thresholds.0 <- .call(
            suggestRnaQcThresholds,
            list(metrics$main),
            list(block=block),
            more.suggest.args
        )
        if (is.null(thresholds)) {
            thresholds <- thresholds.0
        }  else {
            names(thresholds)[names(thresholds) == "subset.proportion"] <- "subsets"
            thresholds <- .substituteThresholds(thresholds.0, thresholds)
        }
    }

    keep <- filterRnaQcMetrics(thresholds, metrics$main, block=block)

    df <- formatComputeRnaQcMetricsResult(metrics$main, flatten=flatten)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::colData(x) <- .cbindOverwrite(SummarizedExperiment::colData(x), df)

    if (!is.null(altexp.proportions)) {
        for (ae.name in names(metrics$altexp)) {
            ae.df <- formatComputeRnaQcMetricsResult(metrics$altexp[[ae.name]], flatten=flatten)
            colnames(ae.df) <- paste0(output.prefix, colnames(ae.df))

            ae.se <- SingleCellExperiment::altExp(x, ae.name)
            ae.cd <- SummarizedExperiment::colData(ae.se) <- .cbindOverwrite(SummarizedExperiment::colData(ae.se), ae.df)
            SingleCellExperiment::altExp(x, ae.name) <- ae.se
        }
    }

    if (!is.null(meta.name)) {
        names(thresholds)[names(thresholds) == "subsets"] <- "subset.proportion"
        metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    if (filter.cells) {
        x <- .delayifyAssays(x)
        x <- x[,keep]
    }

    x
}

.substituteThresholds <- function(oldt, newt, sofar = NULL) {
    for (x in intersect(names(oldt), names(newt))) {
        if (x == "block.ids") {
            next
        }
        oldx <- oldt[[x]]
        newx <- newt[[x]]
        if (is.list(oldx)) {
            oldx <- .substituteThresholds(oldx, newx, sofar = c(sofar, x))
        } else if (is.null(names(newx))) {
            if (length(newx) != 1) {
                stop("expected a scalar at 'thresholds$", paste(sofar, collapse="$"), "'")
            }
            oldx[] <- newx
        } else {
            if (!identical(names(newx), names(oldx))) {
                stop("unexpected names for 'thresholds$", paste(sofar, collapse="$"), "'")
            }
            oldx <- newx
        }
        oldt[[x]] <- oldx
    }
    oldt
}

#' @export
#' @rdname quickRnaQc.se
computeRnaQcMetricsWithAltExps <- function(x, subsets, altexp.proportions, num.threads = 1, assay.type = "counts") {
    metrics <- computeRnaQcMetrics(SummarizedExperiment::assay(x, assay.type), subsets, num.threads=num.threads)

    # Adding more proportions from the alternative experiments, mostly for spike-ins.
    altexp.collected <- list()
    if (!is.null(altexp.proportions)) {
        altexp.proportions <- .sanitizeAltexpAssays(
            altexp.proportions,
            all.altexps=SingleCellExperiment::altExpNames(x),
            default.assay.type=assay.type
        )

        for (ae.name in names(altexp.proportions)) {
            alt.assay <- SummarizedExperiment::assay(SingleCellExperiment::altExp(x, ae.name), altexp.proportions[[ae.name]])
            alt.metrics <- computeRnaQcMetrics(alt.assay, subsets=list(), num.threads=num.threads)
            altexp.collected[[ae.name]] <- alt.metrics
            metrics$subsets[[ae.name]] <- alt.metrics$sum/(metrics$sum + alt.metrics$sum)
        }
    }

    list(main=metrics, altexp=altexp.collected)
}

#' @export
#' @rdname quickRnaQc.se
#' @importFrom S4Vectors DataFrame cbind
formatComputeRnaQcMetricsResult <- function(compute.res, flatten = TRUE) {
    if (flatten) {
        extra <- compute.res$subsets
        colnames(extra) <- sprintf("subset.proportion.%s", colnames(extra))
        cbind(DataFrame(sum=compute.res$sum, detected=compute.res$detected), extra)
    } else {
        colnames(compute.res)[colnames(compute.res) == "subsets"] <- "subset.proportion"
        compute.res
    }
}
