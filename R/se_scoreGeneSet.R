#' Score a gene set in a SummarizedExperiment
#'
#' Compute a gene set activity score for each cell based on the expression values of the genes in the set,
#' by calling \code{\link{scoreGeneSet}} on an assay of a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param set Vector containing the gene set, see \code{?\link{scoreGeneSet}} for details.
#' @param block Block assignment for each cell, passed to \code{\link{scoreGeneSet}}.
#' @param num.threads Number of threads for \code{\link{scoreGeneSet}}.
#' @param more.score.args Named list of further arguments to pass to \code{\link[scrapper]{scoreGeneSet}}.
#' @param assay.type Integer or string specifying the relevant assay in \code{x}, usually containing log-normalized expression values.
#'
#' @return List containing \code{scores}, a numeric vector of the gene set scores across all cells in \code{x};
#' and \code{weights}, a numeric vector of weights for all genes in \code{set}.
#'
#' @author Aaron Lun
#' @examples
#' # Defining a gene set of oligodendrocyte genes.
#' library(org.Mm.eg.db)
#' oligo.set <- select(org.Mm.eg.db, keytype="GO", keys="GO:0048709", columns="SYMBOL")
#' oligo.set <- unique(oligo.set$SYMBOL)
#'
#' sce <- getTestRnaData.se("norm")
#' oligo.scores <- scoreGeneSet.se(sce, oligo.set)
#' summary(oligo.scores$scores)
#'
#' @export
scoreGeneSet.se <- function(
    x,
    set,
    block = NULL,
    num.threads = 1,
    more.score.args = list(),
    assay.type = "logcounts"
) {
    out <- .call(
        scoreGeneSet,
        list(SummarizedExperiment::assay(x, assay.type)),
        list(set=set, block=block, num.threads=num.threads),
        more.score.args
    )

    names(out$scores) <- colnames(x)
    w <- out$weights$weight
    names(w) <- rownames(x)[out$weights$row]
    out$weights <- w

    out
}
