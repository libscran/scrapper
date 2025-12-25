#' Principal components analysis of a Summarizedexperiment
#'
#' Compact and denoise the dataset by performing PCA on the (log-)normalized expression matrix,
#' by calling \code{\link{runPca}} on an assay of a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param features Integer, logical or character vector containing the features of interest to use in the PCA.
#' For RNA data, this is typically the \code{hvg} vector added by \code{\link{chooseRnaHvgs.se}}.
#' If \code{NULL}, all available features are used.
#' @param number Number of PCs to retain, passed to \code{\link{runPca}}.
#' @param num.threads Number of threads for the PCA, passed to \code{\link{runPca}}.
#' @param more.pca.args Named list of additional arguments to pass to \code{\link{runPca}}.
#' @param assay.type Integer or string specifying the assay of \code{x} to be used for PCA.
#' This is typically the log-normalized expression matrix created by \code{\link{normalizeRnaCounts.se}}.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the PC scores.
#' @param meta.name String containing the name of the \code{link[S4Vectors]{metadata}} entry in which to store other PCA statistics.
#' @param delayed.transpose Logical scalar indicating whether to delay the transposition when storing coordinates in the \code{\link[SingleCellExperiment]{reducedDims}}.
#'
#' @return \code{x} is returned with the principal component scores in the \code{reducedDim}.
#' Additional outputs (e.g., rotation matrix, variance explained) are stored in the \code{\link[S4Vectors]{metadata}}.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("hvg")
#' sce <- runPca.se(sce, rowData(sce)$hvg)
#' dim(reducedDim(sce, "PCA"))
#' plot(metadata(sce)$PCA$variance.explained / metadata(sce)$PCA$total.variance)
#'
#' @export
#' @importFrom methods is as
#' @importFrom DelayedArray DelayedArray
runPca.se <- function(
    x,
    features,
    number = 25,
    block = NULL,
    num.threads = 1,
    more.pca.args = list(),
    assay.type = "logcounts",
    output.name = "PCA",
    meta.name = "PCA",
    delayed.transpose = FALSE
) {
    y <- SummarizedExperiment::assay(x, assay.type)
    if (!is.null(features)) {
        y <- DelayedArray(y)[features,,drop=FALSE] # ensure that no copy is made.
    }

    out <- .call(
        scrapper::runPca,
        list(y),
        list(number=number, block=block, num.threads=num.threads),
        more.pca.args
    )

    if (!is(x, "SingleCellExperiment")) {
        loadNamespace("SingleCellExperiment")
        x <- as(x, "SingleCellExperiment")
    }
    x <- .addTransposedReddim(x, output.name, out$components, delayed.transpose)

    if (!is.null(meta.name)) {
        out$components <- NULL
        S4Vectors::metadata(x)[[meta.name]] <- out
    }

    x
}
