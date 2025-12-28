#' Aggregate expression across cells in a SummarizedExperiment
#' 
#' Aggregate expression values across groups of cells for each gene,
#' by calling \code{\link{aggregateAcrossCells}} on an assay in a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param factors List or data frame (or their equivalents from \pkg{S4Vectors}) containing grouping factors, see \code{\link{aggregateAcrossCells}} for more details.
#'
#' Alternatively, an atomic vector or factor representing a single variable.
#' @param num.threads Number of threads, passed to \code{\link{aggregateAcrossCells}}.
#' @param more.aggr.args Named list of additional arguments to pass to \code{\link{aggregateAcrossCells}}.
#' @param assay.type Integer or string specifying the assay of \code{x} to be aggregated.
#' @param output.prefix String specifying a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns storing the factor combinations. 
#' If \code{NULL}, no prefix is added.
#' @param counts.name String specifying the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the cell count for each factor combination. 
#' If \code{NULL}, the cell counts are not reported.
#' @param meta.name String specifying the name of the \code{\link[S4Vectors]{metadata}} entry in which to store additional outputs like the combination indices.
#' If \code{NULL}, additional outputs are not reported.
#' @param include.coldata Logical scalar indicating whether to add the aggregated \code{colData} from \code{x} to the output.
#' @param more.coldata.args Named list of additional arguments to pass to \code{aggregateColData}.
#' Only relevant if \code{include.coldata=TRUE}.
#' @param altexps Unnamed integer or character vector specifying the indices/names of alternative experiments of \code{x} to aggregate.
#' The aggregated assay from each alternative experiment is determined by \code{assay.type}.
#'
#' Alternatively, this may be a named integer or character vector.
#' Each name specifies an alternative experiment while each value is the index/name of the assay to aggregate from that experiment.
#'
#' Only relevant if \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param copy.altexps Logical scalar indicating whether to copy the \code{colData} and \code{metadata} of the output SingleCellExperiment into each of its alternative experiments.
#' @param coldata \link[S4Vectors]{DataFrame} of column data, containing one row for each cell.
#' @param index Integer vector containing the index of the factor combination to which each cell in \code{coldata} was assigned.
#' @param number Integer specifying the total number of unique factor combinations.
#' All elements of \code{index} should be less than \code{number}.
#' @param only.atomic Logical scalar specifying whether to skip non-atomic, non-factor columns.
#' @param placeholder Placeholder value to store in the output column when a factor combination does not have a single unique value. 
#'
#' @return
#' For \code{aggregateAcrossCells.se}, a SummarizedExperiment is returned where each column corresponds to a factor combination.
#' Each row corresponds to a gene in \code{x} and the \code{\link[SummarizedExperiment]{rowData}} is taken from \code{x}.
#' The assays contain the sum of counts (\code{"sums"}) and the number of detected cells (\code{"detected"}) in each combination for each gene.
#' The \code{colData} contains:
#' \itemize{
#' \item The factor combinations, with column names prefixed by \code{output.prefix}.
#' \item The cell count for each combination, named by \code{counts.name}.
#' \item Additional \code{colData} from \code{x} if \code{include.coldata=TRUE}.
#' This is aggregated with \code{aggregateColData} on the combination indices.
#' }
#' The metadata contains a list named as \code{meta.name}, containing a \code{index} integer vector of length equal to the number of cells in \code{x}.
#' Each entry of this vector is an index of the factor combination (i.e., column of the output object) to which the corresponding cell was assigned.
#'
#' If \code{altexps} is specified, a SingleCellExperiment is returned instead.
#' The same aggregation for the main experiment is applied to each alternative experiment.
#' If \code{copy.altexps=TRUE}, the \code{colData} for each alternative experiment will contain a copy of the factor combinations and cell counts,
#' and the \code{metadata} will contain a copy of the index vector.
#'
#' For \code{aggregateColData}, a \link[S4Vectors]{DataFrame} is returned with number of rows equal to \code{number}.
#' Each atomic or factor column in \code{coldata} is represented by a column in the output DataFrame.
#' In each column, the \code{j}-th entry is equal to the unique value of all rows where \code{index == j},
#' or \code{placeholder} if there is not exactly one unique value.
#' If \code{only.atomic=FALSE}, any non-atomic/non-factor columns of \code{coldata} are represented in the output DataFrame by a vector of \code{placeholder} values.
#' If \code{only.atomic=TRUE}, any non-atomic/non-factor columns of \code{coldata} are skipped.
#' 
#' @author Aaron Lun
#' @examples
#' library(SingleCellExperiment)
#' sce <- getTestRnaData.se("start")
#' aggr <- aggregateAcrossCells.se(sce, sce$level1class)
#' head(assay(aggr))
#' colData(aggr)
#'
#' # We can also aggregate within alternative experiments.
#' aggr2 <- aggregateAcrossCells.se(sce, sce$level1class, altexps="ERCC")
#' head(assay(altExp(aggr2, "ERCC")))
#'
#' @export
#' @importFrom methods is
#' @importFrom S4Vectors cbind metadata metadata<-
#' @importClassesFrom S4Vectors List DataFrame
aggregateAcrossCells.se <- function(
    x,
    factors,
    num.threads = 1,
    more.aggr.args = list(),
    assay.type = "counts",
    output.prefix = "factor.",
    counts.name = "counts",
    meta.name = "aggregated",
    include.coldata = TRUE,
    more.coldata.args = list(),
    altexps = NULL,
    copy.altexps = FALSE
) {
    if (is.list(factors) || is.data.frame(factors) || is(factors, "List") || is(factors, "DataFrame")) {
        # this is fine.
    } else {
        factors <- list(factors)
        names(factors) <- make.names(1) # mimic what happens in combineFactors.
    }

    out <- .call(
        aggregateAcrossCells,
        list(SummarizedExperiment::assay(x, assay.type)),
        list(factors=factors, num.threads=num.threads),
        more.aggr.args
    )

    CON <- SummarizedExperiment::SummarizedExperiment
    if (length(altexps)) {
        CON <- SingleCellExperiment::SingleCellExperiment
    }
    se <- CON(out[c("sums", "detected")], rowData=SummarizedExperiment::rowData(x))

    common.cd <- out$combinations
    colnames(common.cd) <- paste0(output.prefix, colnames(common.cd))
    if (!is.null(counts.name)) {
        common.cd[[counts.name]] <- out$counts
    }
    SummarizedExperiment::colData(se) <- common.cd
    if (include.coldata) {
        aggr.cd <- do.call(aggregateColData, c(list(SummarizedExperiment::colData(x), out$index, number=nrow(out$combinations)), more.coldata.args))
        SummarizedExperiment::colData(se) <- cbind(SummarizedExperiment::colData(se), aggr.cd)
    }

    if (!is.null(meta.name)) {
        metadata(se)[[meta.name]] <- list(index=out$index)
    }

    if (length(altexps)) {
        SingleCellExperiment::mainExpName(se) <- SingleCellExperiment::mainExpName(x)
        altexps <- .sanitizeAltexpAssays(altexps, all.altexps=SingleCellExperiment::altExpNames(x), default.assay.type=assay.type)

        for (ae in names(altexps)) {
            ae.se <- Recall(
                SingleCellExperiment::altExp(x, ae),
                out$index,
                num.threads=num.threads,
                more.aggr.args=more.aggr.args,
                assay.type=altexps[[ae]],
                altexps=NULL,
                output.prefix=NULL,
                counts.name=NULL,
                meta.name=NULL,
                include.coldata=include.coldata
            )

            ae.cd <- SummarizedExperiment::colData(ae.se)[,-1,drop=FALSE] # remove uninteresting factor combination
            if (copy.altexps) {
                ae.cd <- cbind(common.cd, ae.cd)
            }
            SummarizedExperiment::colData(ae.se) <- ae.cd

            if (copy.altexps) {
                metadata(ae.se) <- metadata(se)
            }
            SingleCellExperiment::altExp(se, ae) <- ae.se
        }
    }

    se
}

#' @export
#' @rdname aggregateAcrossCells.se
#' @importFrom S4Vectors make_zero_col_DFrame
aggregateColData <- function(coldata, index, number, only.atomic = TRUE, placeholder = NA) {
    collected <- make_zero_col_DFrame(nrow=number)
    index <- factor(index, seq_len(number))

    for (cn in colnames(coldata)) {
        curcol <- coldata[[cn]]
        if (!is.atomic(curcol) && !is.factor(curcol)) {
            if (!only.atomic) {
                collected[[cn]] <- rep(placeholder, number)
            }
            next
        }

        grouped <- split(curcol, index)
        alloc <- rep(curcol[1], number)
        for (i in seq_along(grouped)) {
            u <- unique(grouped[[i]])
            if (length(u) != 1L) {
                alloc[i] <- placeholder
            } else {
                alloc[i] <- u
            }
        }

        collected[[cn]] <- alloc
    }

    collected
}
