#' Score marker genes in a SummarizedExperiment
#'
#' Identify candidate marker genes based on effect sizes from pairwise comparisons between groups of cells,
#' by calling \code{\link{scoreMarkers}} on an assay of a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param groups Group assignment for each cell, passed to \code{\link{scoreMarkers}}.
#' @param block Block assignment for each cell, passed to \code{\link{scoreMarkers}}.
#' @param num.threads Number of threads for marker scoring, passed to \code{\link{scoreMarkers}}.
#' @param more.marker.args Named list of additional arguments to pass to \code{\link{scoreMarkers}}.
#' @param assay.type Integer or string specifying the assay to use for differential comparisons, usually containing log-normalized expression values.
#' @param extra.columns \link[S4Vectors]{DataFrame} containing extra columns to add each DataFrame.
#' This should have the same number of rows as \code{x}.
#' For \code{scoreMarkers.se}, this may also be a character vector specifying the columns of \code{\link[SummarizedExperiment]{rowData}} to be added.
#' @param order.by String specifying the column to order each DataFrame by.
#' Alternatively \code{TRUE}, a column is automatically chosen from the effect size summaries.
#' If \code{NULL} or \code{FALSE}, no ordering is performed.
#' @param marker.res List containing the result of \code{\link{scoreMarkers}}.
#' @param marker.df DataFrame containing the marker statistics for a single group.
#' @param columns Character vector of the columns to retain in the preview.
#' This may be named, in which the names are used as the column names.
#' @param include.order.by Boolean indicating whether the column specified by \code{order.by} should be included in the output DataFrame.
#' @param rows Integer specifying the number of rows to show.
#' If \code{NULL}, all rows are returned.
#' 
#' @return 
#' For \code{scoreMarkers.se} and \code{formatScoreMarkersResult}, a \link[S4Vectors]{List} of \link[S4Vectors]{DataFrame}s is returned.
#' Each DataFrame corresponds to a unique group in \code{groups}.
#' Each row contains statistics for a gene in \code{x}, with the following columns:
#' \itemize{
#' \item \code{mean}, the mean expression in the current group.
#' \item \code{detected}, the proportion of cells with detected expression in the current group.
#' \item \code{<effect>.<summary>}, a summary statistic for an effect size,
#' e.g., \code{cohens.d.mean} contains the mean Cohen's d across comparisons involving the current group.
#' }
#'
#' For \code{previewMarkers}, a DataFrame is returned containing the specified columns and rows.
#'
#' @author Aaron Lun
#'
#' @examples
#' sce <- getTestRnaData.se("cluster")
#' markers <- scoreMarkers.se(sce, sce$clusters)
#' previewMarkers(markers[["1"]], c(effect="cohens.d.mean"))
#'
#' @export
scoreMarkers.se <- function(
    x,
    groups,
    block = NULL,
    num.threads = 1,
    more.marker.args = list(),
    assay.type = "logcounts",
    extra.columns = NULL,
    order.by = TRUE
) {
    res <- .call(
        scoreMarkers,
        list(SummarizedExperiment::assay(x, assay.type), groups=groups),
        list(block=block, num.threads=num.threads),
        more.marker.args
    )

    if (is.character(extra.columns)) {
        extra.columns <- SummarizedExperiment::rowData(x)[,extra.columns,drop=FALSE]
    }
    formatScoreMarkersResult(res, extra.columns=extra.columns, order.by=order.by)
}

.guessDimnames <- function(marker.res) {
    for (n in names(marker.res)) {
        current <- marker.res[[n]]
        if (is.matrix(current)) {
            return(list(nrow=nrow(current), rownames=rownames(current), groups=colnames(current)))
        } else if (is.data.frame(current)) {
            return(list(nrow=nrow(current), rownames=rownames(current), groups=NULL))
        } else if (is.list(current)) {
            out <- .guessDimnames(current)
            if (!is.null(out)) {
                out$groups <- names(current)
                return(out)
            }
        } else {
            stop("unknown type '", typeof(current), "'")
        }
    }
    return(NULL)
}

.findOrderBy <- function(marker.df, order.by) {
    if (isTRUE(order.by)) {
        # Find something decent to use for ordering.
        for (summ in c("mean", "median", "min.rank", "min", "max")) {
            for (eff in c("cohens.d", "auc", "delta.mean", "delta.detected")) {
                candidate <- paste0(eff, ".", summ)
                if (candidate %in% names(marker.df)) {
                    return(candidate)
                }
            }
        }
        return(NULL)
    } else if (isFALSE(order.by)) {
        return(NULL)
    } else {
        # No-op if it was already NULL or a string.
        return(order.by)
    }
}

#' @export
#' @rdname scoreMarkers.se
#' @importFrom S4Vectors List make_zero_col_DFrame cbind
formatScoreMarkersResult <- function(marker.res, extra.columns = NULL, order.by = TRUE) {
    effect.sizes <- c("cohens.d", "auc", "delta.mean", "delta.detected")
    summaries <- c("min", "mean", "median", "max", "quantile", "min.rank")

    out <- .guessDimnames(marker.res)
    if (is.null(out)) {
        stop("could not determine dimnames from 'marker.res'")
    }

    output <- List()
    check.order.by <- FALSE

    for (group in out$groups) {
        current <- make_zero_col_DFrame(out$nrow)
        rownames(current) <- out$rownames
        if (!is.null(extra.columns)) {
            current <- cbind(current, extra.columns)
        }

        if ("mean" %in% names(marker.res)) {
            current$mean <- marker.res$mean[,group]
        }
        if ("detected" %in% names(marker.res)) {
            current$detected <- marker.res$detected[,group]
        }

        for (eff in effect.sizes) {
            eff.all <- marker.res[[eff]]
            if (is.null(eff.all)) {
                next
            }

            eff.df <- eff.all[[group]]
            for (summ in summaries) {
                eff.summ <- eff.df[[summ]]
                if (is.null(eff.summ)) {
                    next
                }

                if (is.data.frame(eff.summ)) {
                    colnames(eff.summ) <- sprintf("%s.%s.%s", eff, summ, colnames(eff.summ))
                    current <- cbind(current, eff.summ)
                } else {
                    current[[paste0(eff, ".", summ)]] <- eff.summ
                }
            }
        }

        if (!check.order.by) {
            order.by <- .findOrderBy(current, order.by)
            check.order.by <- TRUE
        }
        if (!is.null(order.by)) {
            dec <- !endsWith(order.by, ".min.rank")
            current <- current[order(current[[order.by]], decreasing=dec),,drop=FALSE]
        }

        output[[group]] <- current
    }

    output
}

#' @export
#' @rdname scoreMarkers.se
#' @importFrom utils head
previewMarkers <- function(marker.df, columns = c("mean", "detected", lfc="delta.mean.mean"), rows = 10, order.by = NULL, include.order.by = !is.null(order.by)) {
    if (is.null(columns)) {
        columns <- as.character(0)
    }

    order.by <- .findOrderBy(marker.df, order.by)
    if (include.order.by) {
        columns <- c(columns, order.by)
    }

    columns <- columns[columns %in% colnames(marker.df)]
    df <- marker.df[,columns,drop=FALSE]
    if (!is.null(names(columns))) {
        replace <- names(columns) != ""
        colnames(df)[replace] <- names(columns)[replace]
    }

    if (!is.null(order.by)) {
        dec <- !endsWith(order.by, ".min.rank")
        o <- order(marker.df[[order.by]], decreasing=dec) # use marker.df in case the ordering statistic was renamed.
        if (!is.null(rows)) {
            o <- head(o, rows)
        }
        df <- df[o,,drop=FALSE]
    } else {
        if (!is.null(rows)) {
            df <- head(df, rows)
        }
    }

    df
}
