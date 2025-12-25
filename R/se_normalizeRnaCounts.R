#' Normalize RNA counts in a SummarizedExperiment
#'
#' Compute (log-)normalized expression values after performing scaling normalization of an RNA count matrix.
#' This calls \code{\link{normalizeCounts}} on an assay of a \link[SummarizedExperiment]{SummarizedExperiment},
#' after centering the size factors with \code{\link{centerSizeFactors}}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param size.factors Numeric vector of length equal to the number of columns of \code{x},
#' containing the size factor for each cell in \code{x}.
#' If \code{NULL}, this defaults to the column sums of the count matrix in \code{x}.
#' @param center Logical scalar indicating whether to center the \code{size.factors},
#' see \code{?\link[scrapper]{centerSizeFactors}} for more details.
#' @param block Block assignments for each cell, passed to \code{\link[scrapper]{centerSizeFactors}}.
#' @param mode How to center size factors in different blocks, see \code{?\link[scrapper]{centerSizeFactors}} for more details.
#' @param log Whether to log-transform the normalized expression values, see \code{?\link[scrapper]{normalizeCounts}} for more details.
#' @param pseudo.count The pseudo-count for log-transformation, see \code{?\link[scrapper]{normalizeCounts}} for more details.
#' @param assay.type Integer or string specifying the assay of \code{x} with the count matrix.
#' @param output.name String containing the name of the assay to store the normalized matrix.
#' @param factor.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the size factors in the output object.
#' If \code{NULL}, the size factors are not stored. 
#'
#' @return \code{x} is returned with a new assay containing the (log-)normalized matrix.
#' Size factors are also stored in the \code{\link[SummarizedExperiment]{colData}}.
#'
#' @author Aaron Lun
#'
#' @examples
#' sce <- getTestRnaData.se("qc")
#' sce <- normalizeRnaCounts.se(sce, size.factors=sce$sum)
#' assayNames(sce)
#' summary(sizeFactors(sce))
#'
#' @export
#' @importFrom Matrix colSums
normalizeRnaCounts.se <- function(
    x,
    size.factors = NULL,
    center = TRUE,
    block = NULL,
    mode = "lowest",
    log = TRUE,
    pseudo.count = 1,
    assay.type = "counts",
    output.name = "logcounts",
    factor.name = "sizeFactor"
) {
    y <- SummarizedExperiment::assay(x, assay.type)

    if (is.null(size.factors)) {
        size.factors <- colSums(y)
    }
    if (center) {
        size.factors <- centerSizeFactors(size.factors, block=block, mode=mode)
    }

    SummarizedExperiment::assay(x, output.name) <- normalizeCounts(y, size.factors=size.factors, log=log, pseudo.count=pseudo.count)
    if (!is.null(factor.name)) {
        SummarizedExperiment::colData(x)[[factor.name]] <- size.factors
    }

    x
}
