#' Aggregate expression across gene sets in a SummarizedExperiment
#'
#' Aggregate expression values across sets of genes for each cell,
#' by calling \code{\link{aggregateAcrossGenes}} on an assay in a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param sets List of gene sets, see \code{\link{aggregateAcrossGenes}} for more details.
#'
#' Alternatively, \code{sets} may be a \link[S4Vectors]{List} subclass,
#' in which case the \code{\link[S4Vectors]{mcols}} are used as the \code{\link[SummarizedExperiment]{rowData}} of the output object.
#' Weighted gene sets can be represented by a list of \pkg{DataFrames} or a \link[IRanges]{DataFrameList},
#' where each DataFrame contains two columns, i.e., the gene identities and the associated weights.
#' @param num.threads Number of threads, passed to \code{\link{aggregateAcrossGenes}}.
#' @param more.aggr.args Named list of additional arguments to pass to \code{\link{aggregateAcrossGenes}}.
#' @param assay.type Integer or string specifying the assay of \code{x} to be aggregated.
#' @param output.name String specifying the assay name in the output object.
#' Defaults to \code{assay.type} if it is a string, otherwise \code{"aggregated"}.
#'
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with number of rows equal to the number of gene sets.
#' The lone assay contains the aggregated values for each gene set for all cells.
#' The \code{\link[SummarizedExperiment]{colData}} is the same as that of \code{x}.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("norm")
#'
#' library(org.Mm.eg.db)
#' set.df <- select(
#'     org.Mm.eg.db,
#'     keytype="GO",
#'     keys=c(
#'         "GO:0048709", # oligodendrocyte differentiation 
#'         "GO:0048699", # neuron development
#'         "GO:0048143"  # astrocyte activation 
#'     ),
#'     columns="SYMBOL"
#' )
#' sets <- splitAsList(set.df$SYMBOL, set.df$GO)
#'
#' aggregated <- aggregateAcrossGenes.se(sce, sets)
#' aggregated
#' assay(aggregated)[,1:10]
#'
#' @export
#' @importFrom methods is
#' @importFrom S4Vectors mcols
aggregateAcrossGenes.se <- function(
    x,
    sets,
    num.threads = 1,
    more.aggr.args = list(),
    assay.type = "logcounts",
    output.name = NULL
) {
    rd <- NULL
    if (is(sets, "List")) {
        rd <- mcols(sets)
        if (!is.null(rd)) {
            rownames(rd) <- names(sets) # just in case
        }
        sets <- as.list(sets)
    }

    for (i in seq_along(sets)) {
        if (is(sets[[i]], "List")) {
            sets[[i]] <- as.list(sets[[i]])
        }
    }

    vecs <- .call(
        aggregateAcrossGenes,
        list(SummarizedExperiment::assay(x, assay.type)),
        list(sets=sets, num.threads=num.threads),
        more.aggr.args
    )
    if (length(vecs)) {
        mat <- do.call(rbind, vecs)
    } else {
        mat <- matrix(0, 0, ncol(x))
    }

    assays <- list(mat)
    if (!is.null(output.name)) {
        names(assays) <- output.name
    } else if (is.character(assay.type)) {
        names(assays) <- assay.type
    } else {
        names(assays) <- "aggregated"
    }

    out <- SummarizedExperiment::SummarizedExperiment(assays, colData=SummarizedExperiment::colData(x))
    if (!is.null(rd)) {
        SummarizedExperiment::rowData(out) <- rd
    }
    out
}
