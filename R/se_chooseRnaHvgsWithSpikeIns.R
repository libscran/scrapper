#' Choose highly variable genes based on spike-ins
#'
#' Fit a mean-variance trend to spike-in transcripts, and use this to estimate the technical noise for endogenous genes of similar abundance.
#' Highly variable genes are then selected based on the biological component of variation.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} or one of its subclasses.
#' The main experiment should contain count data for endogenous genes, where rows correspond to genes and columns correspond to cells.
#' @param spike.altexp Integer or string specifying the name/index of the alternative experiment containing the spike-in data.
#' The assay to use is determined by \code{assay.type}.
#'
#' Alternatively, \code{spike.altexp} may be a named integer or character vector of length 1.
#' The name specifies an alternative experiment while the value is the index/name of the assay to use from that experiment.
#' @param more.endogenous.var.args Named list of additional arguments to pass to \code{\link{modelGeneVariances}} for endogenous genes.
#' Note that all trend-fitting arguments will be ignored as the mean-variance trend is obtained from the spike-in transcripts.
#' @param more.spike.var.args Named list of additional arguments to pass to \code{\link{modelGeneVariances}} for spike-in transcripts.
#' @inheritParams chooseRnaHvgs.se
#'  
#' @return
#' \code{x} is returned with the per-gene variance modelling statistics for endogenous genes added to its \code{rowData}.
#' This is equivalent to the columns returned by \code{\link{modelGeneVariances}}, except that:
#' \itemize{
#' \item \code{fitted} now contains the interpolated values from the mean-variance trend fitted to the spike-in transcripts.
#' This represents the technical component of variation for each gene.
#' \item \code{residuals} now contains the difference of each gene's variance from the trend.
#' This represents the biological component of variation for each gene.
#' }
#' The \code{hvg} column in the \code{rowData} indicates whether a gene was chosen as a HVG based on the largest \code{residuals}.
#' If \code{include.per.block=TRUE} and \code{block} is specified, the per-block statistics are stored as a nested DataFrame in the \code{per.block} column.
#'
#' Per-transcript modelling statistics are also added to the \code{rowData} of the alternative experiment containing the spike-in data.
#' This has all of the same columns returned by \code{\link{modelGeneVariances}}.
#'
#' @author Aaron Lun
#'
#' @details
#' It is generally expected that normalization has been performed with \code{\link{normalizeRnaCountsWithSpikeIns.se}}.
#' This ensures that the means are comparable between endogenous genes and spike-in transcripts.
#'
#' @seealso
#' \code{\link{chooseRnaHvgs.se}}, for a simpler function when spike-ins are not available.
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("qc")
#' sce <- normalizeRnaCountsWithSpikeIns.se(sce, "ERCC")
#'
#' sce <- chooseRnaHvgsWithSpikeIns.se(sce, "ERCC")
#' summary(rowData(sce)$hvg)
#'
#' plot(rowData(sce)$means, rowData(sce)$variances, col=factor(rowData(sce)$hvg))
#' spike.rd <- rowData(altExp(sce, "ERCC"))
#' points(spike.rd$means, spike.rd$variances, col="dodgerblue", pch=4)
#' curve(approxfun(spike.rd$means, spike.rd$fitted)(x), col="dodgerblue", add=TRUE)
#' 
#' @export
chooseRnaHvgsWithSpikeIns.se <- function(
    x,
    spike.altexp,
    block = NULL,
    num.threads = 1,
    more.endogenous.var.args = list(),
    more.spike.var.args = list(),
    top = 4000,
    more.choose.args = list(),
    assay.type = "logcounts",
    output.prefix = NULL,
    include.per.block = FALSE
) {
    spike.altexp <- .sanitizeAltexpAssays(
        spike.altexp,
        all.altexps=SingleCellExperiment::altExpNames(x),
        default.assay.type=assay.type
    )

    common.var.args <- list(block=block, num.threads=num.threads)

    endog.var <- .call(
        modelGeneVariances,
        list(SummarizedExperiment::assay(x, assay.type), fit.trend=FALSE),
        common.var.args,
        more.endogenous.var.args
    )

    spike.x <- SingleCellExperiment::altExp(x, names(spike.altexp)[1])
    spike.var <- .call(
        modelGeneVariances,
        list(SummarizedExperiment::assay(spike.x, spike.altexp[[1]])),
        common.var.args,
        more.endogenous.var.args
    )

    # Adding zero to force a linear to zero on the left edge.
    endog.var$statistics$fitted <- approx(c(spike.var$statistics$means, 0), c(spike.var$statistics$fitted, 0), xout=endog.var$statistics$means, rule=2)$y
    endog.var$statistics$residuals <- endog.var$statistics$variances - endog.var$statistics$fitted

    hvg.index <- .call(
        chooseHighlyVariableGenes,
        list(endog.var$statistics$residuals),
        list(top=top, larger=TRUE),
        more.choose.args
    )

    endog.df <- formatModelGeneVariancesResult(endog.var, choose.res=hvg.index, include.per.block=include.per.block)
    colnames(endog.df) <- paste0(output.prefix, colnames(endog.df))
    SummarizedExperiment::rowData(x) <- cbind(SummarizedExperiment::rowData(x), endog.df)

    spike.df <- formatModelGeneVariancesResult(spike.var, choose.res=NULL, include.per.block=include.per.block)
    colnames(spike.df) <- paste0(output.prefix, colnames(spike.df))
    SummarizedExperiment::rowData(spike.x) <- cbind(SummarizedExperiment::rowData(spike.x), spike.df)
    SingleCellExperiment::altExp(x, names(spike.altexp)[1]) <- spike.x

    x
}
