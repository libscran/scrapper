#' Normalize RNA and spike-in counts
#'
#' Compute (log-)normalized expression values for endogenous genes and spike-in transcripts after scaling normalization.
#' All sets of size factors are centered with \code{\link{centerSpikeInFactors}} prior to calling \code{\link{normalizeCounts}} on the relevant assays.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} or one of its subclasses.
#' The main experiment should contain count data for endogenous genes, where rows correspond to genes and columns correspond to cells.
#' There should also be one or more alternative experiments containing spike-in data, see \code{spike.altexps}.
#' @param spike.altexps Alternative experiments containing spike-in data.
#' This should be an unnamed integer or character vector containing the names/indices of the alternative experiments of interest.
#' The assay to use from each alternative experiment is determined by \code{assay.type}.
#'
#' Alternatively, \code{spike.altexps} may be a named integer or character vector.
#' Each name specifies an alternative experiment while each value is the index/name of the assay to use from that experiment.
#' @param endogenous.factors Numeric vector of length equal to \code{ncol(x)}, containing the size factors for the endogenous genes.
#' If \code{NULL}, this defaults to the column sums of the relevant count matrix in \code{x}.
#' Ignored if \code{use.spike.ins.for.endogenous.factors} is not \code{FALSE}.
#' @param spike.factors Named list of numeric vectors containing the size factors for each spike-in set.
#' Each entry corresponds to a spike-in set ad be named after its alternative experiment in \code{spike.altexps}.
#' Each vector should be of length equal to \code{ncol(x)}. 
#' If \code{NULL} or a spike-in set is not present in \code{spike.factors}, the column sums of the relevant count matrix is used instead.
#' @inheritParams normalizeRnaCounts.se
#' @param num.threads Integer specifying the number of threads for computing column sums.
#' Only used if precomputed factors are not available in \code{endogenous.factors} or \code{spike.factors}.
#' @param use.spike.ins.for.endogenous Boolean indicating whether to use spike-in factors for normalization of endogenous genes.
#' By default, each set of spike-in factors is only used on its corresponding count matrix, while the endogenous factors are used on the endogenous count matrix.
#'
#' If this is set to \code{TRUE}, the size factors from the first spike-in set in \code{spike.altexps} are used to normalize the endogenous counts.
#' This is useful for removing technical differences in scaling while preserving differences in total RNA content between cells.
#'
#' Alternatively, this may be a string specifying the name of the spike-in set from which to obtain size factors for normalization of the endogenous counts.
#' This should refer to one of the alternative experiments in \code{spike.altexps}.
#'
#' @return \code{x} is returned with new assays containing (log-)normalized matrices for endogenous genes and spike-in transcripts.
#' Size factors are also stored in the \code{\link[SummarizedExperiment]{colData}} for each experiment.
#'
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("qc")
#' sce <- normalizeRnaCountsWithSpikeIns.se(sce, "ERCC")
#'
#' summary(sce$sizeFactor)
#' assayNames(sce)
#'
#' summary(altExp(sce, "ERCC")$sizeFactor)
#' assayNames(altExp(sce, "ERCC"))
#'
#' @seealso
#' \code{\link{normalizeRnaCounts.se}}, for simpler normalization when spike-ins are not present.
#' 
#' @export
#' @importFrom beachmat initializeCpp tatami.column.sums
normalizeRnaCountsWithSpikeIns.se <- function(
    x,
    spike.altexps,
    endogenous.factors = NULL,
    spike.factors = NULL,
    use.spike.ins.for.endogenous = FALSE,
    block = NULL,
    mode = "lowest",
    log = TRUE,
    pseudo.count = 1,
    more.norm.args = list(),
    assay.type = "counts",
    output.name = "logcounts",
    factor.name = "sizeFactor",
    num.threads = 1
) {
    spike.altexps <- .sanitizeAltexpAssays(
        spike.altexps,
        all.altexps=SingleCellExperiment::altExpNames(x),
        default.assay.type=assay.type
    )

    spike.xs <- list()
    spike.ys <- list()
    for (sp in names(spike.altexps)) {
        spike.x <- SingleCellExperiment::altExp(x, sp)
        spike.y <- SummarizedExperiment::assay(spike.x, spike.altexps[[sp]])
        if (!(sp %in% names(spike.factors))) {
            spike.factors[[sp]] <- tatami.column.sums(initializeCpp(spike.y), num.threads=num.threads)
        }
        spike.ys[[sp]] <- spike.y
        spike.xs[[sp]] <- spike.x
    }

    endog.y <- SummarizedExperiment::assay(x, assay.type)

    center.args <- list(block=block, mode=mode)
    if (isFALSE(use.spike.ins.for.endogenous)) {
        if (is.null(endogenous.factors)) {
            endogenous.factors <- tatami.column.sums(initializeCpp(endog.y), num.threads=num.threads)
        }
        computed <- do.call(centerSpikeInFactors, c(list(endogenous.factors, spike.factors), center.args))

    } else {
        if (isTRUE(use.spike.ins.for.endogenous)) {
            use.spike.ins.for.endogenous <- names(spike.altexps)[1]
        }
        endogenous.factors <- spike.factors[[use.spike.ins.for.endogenous]] 
        spike.factors[[use.spike.ins.for.endogenous]] <- NULL # avoid unnecessary redundant calculation. 
        computed <- do.call(centerSpikeInFactors, c(list(endogenous.factors, spike.factors), center.args))
        computed$spike.ins[[use.spike.ins.for.endogenous]] <- computed$endogenous
    }

    SummarizedExperiment::assay(x, output.name) <- .call(
        normalizeCounts,
        list(endog.y, size.factors=computed$endogenous),
        list(log=log, pseudo.count=pseudo.count),
        more.norm.args
    )
    if (!is.null(factor.name)) {
        SummarizedExperiment::colData(x)[[factor.name]] <- computed$endogenous 
    }

    for (sp in names(spike.altexps)) {
        spike.x <- spike.xs[[sp]]
        spike.sf <- computed$spike.ins[[sp]]
        SummarizedExperiment::assay(spike.x, output.name) <- .call(
            normalizeCounts,
            list(spike.ys[[sp]], size.factors=spike.sf),
            list(log=log, pseudo.count=pseudo.count),
            more.norm.args
        )
        if (!is.null(factor.name)) {
            SummarizedExperiment::colData(spike.x)[[factor.name]] <- spike.sf
        }
        SingleCellExperiment::altExp(x, sp) <- spike.x
    }

    x
}
