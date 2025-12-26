#' Choose highly variable genes from a SummarizedExperiment 
#'
#' Model the mean-variance relationship across genes and choose highly variable genes (HVGs) based on the residuals of the fitted trend.
#' This calls \code{\link{modelGeneVariances}} on an assay of a \link[SummarizedExperiment]{SummarizedExperiment},
#' and then calls \code{\link{chooseHighlyVariableGenes}} on the residuals.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param block Block assignment for each cell, to pass to \code{\link{modelGeneVariances}}.
#' @param num.threads Number of threads, to pass to \code{\link{modelGeneVariances}}.
#' @param more.var.args Named list of arguments to pass to \code{\link{modelGeneVariances}}.
#' @param top Number of HVGs to choose, to pass to \code{\link{chooseHighlyVariableGenes}}.
#' @param more.choose.args Named list of arguments to pass to \code{\link{chooseHighlyVariableGenes}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the log-normalized expression matrix for the RNA data.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{rowData}} columns containing the output statistics.
#' @param include.per.block Logical scalar indicating whether the per-block statistics should be stored in the output \code{rowData}.
#' Only relevant if \code{block} is specified.
#' @param model.res List returned by \code{\link{modelGeneVariances}}.
#' @param choose.res Integer vector returned by \code{\link{chooseHighlyVariableGenes}}.
#' This may be \code{NULL}, in which case the identities of the HVGs will not be stored.
#'
#' @return
#' For \code{chooseRnaHvgs.se}, \code{x} is returned with the per-gene variance modelling statistics added to its \code{rowData}.
#' The \code{hvg} column in the \code{rowData} indicates whether a gene was chosen as a HVG.
#' If \code{include.per.block=TRUE} and \code{block} is specified, the per-block statistics are stored as a nested DataFrame in the \code{per.block} column.
#'
#' For \code{formatModelGeneVariancesResult}, a \link[S4Vectors]{DataFrame} is returned with the per-gene variance modelling statistics.
#' If \code{choose.res} is provided, a \code{hvg} column is also stored that indicates whether a gene was chosen as a HVG.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("norm")
#' sce <- chooseRnaHvgs.se(sce, more.var.args=list(use.min.width=TRUE))
#' summary(rowData(sce)$hvg)
#'
#' plot(rowData(sce)$means, rowData(sce)$variances, col=factor(rowData(sce)$hvg))
#' curve(approxfun(rowData(sce)$means, rowData(sce)$fitted)(x), col="dodgerblue", add=TRUE)
#' 
#' @export
chooseRnaHvgs.se <- function(
    x, 
    block = NULL,
    num.threads = 1,
    more.var.args = list(),
    top = 4000,
    more.choose.args = list(),
    assay.type = "logcounts",
    output.prefix = NULL,
    include.per.block = FALSE
) {
    info <- .call(
        modelGeneVariances,
        list(SummarizedExperiment::assay(x, assay.type)),
        list(block=block, num.threads=num.threads),
        more.var.args
    )

    hvg.index <- .call(
        chooseHighlyVariableGenes,
        list(info$statistics$residuals),
        list(top=top, larger=TRUE),
        more.choose.args
    )

    df <- formatModelGeneVariancesResult(info, choose.res=hvg.index, include.per.block=include.per.block)
    colnames(df) <- paste0(output.prefix, colnames(df))
    SummarizedExperiment::rowData(x) <- S4Vectors::cbind(SummarizedExperiment::rowData(x), df)
    x
}

#' @export
#' @rdname chooseRnaHvgs.se
formatModelGeneVariancesResult <- function(model.res, choose.res = NULL, include.per.block = FALSE) {
    df <- S4Vectors::DataFrame(model.res$statistics)

    if (include.per.block && !is.null(model.res$per.block)) {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(df))
        for (n in names(model.res$per.block)) {
            tmp[[n]] <- S4Vectors::DataFrame(model.res$per.block[[n]])
        }
        df[["per.block"]] <- tmp
    }

    if (!is.null(choose.res)) {
        is.hvg <- logical(nrow(df))
        is.hvg[choose.res] <- TRUE
        df[["hvg"]] <- is.hvg
    }

    df
}
