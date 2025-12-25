#' Normalize ADT counts in a SummarizedExperiment
#'
#' Compute (log-)normalized expression values after performing scaling normalization of an ADT count matrix.
#' This calls \code{\link{computeClrm1Factors}} on an assay of a \link[SingleCellExperiment]{SingleCellExperiment},
#' centering the subsequent size factors with \code{\link{centerSizeFactors}},
#' and then computing normalized log-expression values with \code{\link{normalizeCounts}}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to antibody-derived tags (ADTs) and columns correspond to cells.
#' @param size.factors Numeric vector of length equal to the number of columns of \code{x},
#' containing the size factor for each cell in \code{x}.
#' If \code{NULL}, this defaults to the output of \code{\link[scrapper]{computeClrm1Factors}}.
#' @param num.threads Number of threads, passed to \code{\link[scrapper]{computeClrm1Factors}}.
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
#' library(SummarizedExperiment)
#' sce <- altExp(getTestAdtData.se("qc"), "ADT")
#' sce <- normalizeAdtCounts.se(sce)
#' assayNames(sce)
#' summary(sizeFactors(sce))
#'
#' @export
normalizeAdtCounts.se <- function(
    x,
    size.factors = NULL,
    num.threads = 1,
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
        size.factors <- computeClrm1Factors(y, num.threads=num.threads)
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
