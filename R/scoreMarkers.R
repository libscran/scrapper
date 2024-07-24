#' Score marker genes
#'
#' Score marker genes for each group using a variety of effect sizes from pairwise comparisons between groups.
#' This includes Cohen's d, the area under the curve (AUC), the difference in the means (delta-mean) and the difference in the proportion of detected cells (delta-detected).
#'
#' @inheritParams modelGeneVariances
#' @param groups A vector specifying the group assignment for each cell in \code{x}.
#' @param threshold Non-negative numeric scalar specifying the minimum threshold on the differences in means (i.e., the log-fold change, if \code{x} contains log-expression values). 
#' This is incorporated into the effect sizes for Cohen's d and the AUC.
#' @param compute.auc Logical scalar indicating whether to compute the AUC.
#' Setting this to \code{FALSE} can improve speed and memory efficiency.
#' @param all.pairwise Logical scalar indicating whether to report the full effects for every pairwise comparison between groups.
#'
#' @return If \code{all.pairwise=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{mean}, a numeric matrix containing the mean expression for each group.
#' Each row is a gene and each column is a group.
#' \item \code{detected}, a numeric matrix containing the proportion of detected cells in each group.
#' Each row is a gene and each column is a group.
#' \item \code{cohens.d}, a list of data frames where each data frame corresponds to a group.
#' Each row of each data frame represents a gene, while each column contains a summary of Cohen's d from pairwise comparisons to all other groups.
#' This includes the \code{min}, \code{mean}, \code{median}, \code{max} and \code{min.rank}.
#' \item \code{auc}, a list like \code{cohens.d} but containing the summaries of the AUCs from each pairwise comparison.
#' Omitted if \code{compute.auc=FALSE}.
#' \item \code{delta.mean}, a list like \code{cohens.d} but containing the summaries of the delta-mean from each pairwise comparison.
#' \item \code{delta.mean}, a list like \code{cohens.d} but containing the summaries of the delta-detected from each pairwise comparison.
#' }
#'
#' If \code{all.pairwise=TRUE}, a list is returned containing:
#' \itemize{
#' \item \code{mean}, a numeric matrix containing the mean expression for each group.
#' Each row is a gene and each column is a group.
#' \item \code{detected}, a numeric matrix containing the proportion of detected cells in each group.
#' Each row is a gene and each column is a group.
#' \item \code{cohens.d}, a 3-dimensional numeric array containing the Cohen's from each pairwise comparison between groups.
#' The first dimension represents the first group, the second dimension represents the second group, and the final dimension represents the gene;
#' the entry \code{[i, j, k]} represents Cohen's d for \code{i} minus \code{j} of gene \code{k}.
#' \item \code{auc}, an array like \code{cohens.d} but containing the AUCs from each pairwise comparison.
#' Omitted if \code{compute.auc=FALSE}.
#' \item \code{delta.mean}, an array like \code{cohens.d} but containing the delta-mean from each pairwise comparison.
#' \item \code{delta.detected}, an array like \code{cohens.d} but containing the delta-detected from each pairwise comparison.
#' }
#'
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' normed <- normalizeCounts(x, size.factors=centerSizeFactors(colSums(x)))
#'
#' g <- sample(letters[1:4], ncol(x), replace=TRUE)
#' scores <- scoreMarkers(normed, g)
#' names(scores)
#' head(scores$mean)
#' head(scores$cohens.d[["a"]])
#'
#' @seealso
#' \url{https://libscran.github.io/scran_markers}, in particular 
#' the \code{score_markers_summary} function (for \code{all.pairwise=FALSE}),
#' the \code{score_markers_pairwise} function (for \code{all.pairwise=TRUE}),
#' and their blocked equivalents \code{score_markers_summary_blocked} and \code{score_markers_pairwise_blocked} (when \code{block} is not \code{NULL}).
#'
#' \code{\link{summarizeEffects}}, to summarize the pairwise effects returned when \code{all.pairwise=TRUE}.
#' @export
#' @importFrom beachmat initializeCpp
scoreMarkers <- function(
    x, 
    groups, 
    block=NULL, 
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    compute.auc=TRUE,
    threshold=0, 
    all.pairwise=FALSE, 
    num.threads=1)
{
    x <- initializeCpp(x)
    groups <- .transformFactor(groups)
    block <- .transformFactor(block)

    args <- list(
        groups=groups$index,
        num_groups=length(groups$names),
        block=block$index,
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        threshold=threshold,
        compute_auc=compute.auc,
        num_threads=num.threads
    )

    if (all.pairwise) {
        output <- do.call(score_markers_pairwise, c(list(x), args))
        dimnames(output$mean) <- dimnames(output$detected) <- list(rownames(x), groups$names)

        for (nm in c("cohens.d", "delta.mean", "delta.detected", "auc")) {
            current <- output[[nm]]
            if (length(current)) {
                dimnames(current) <- list(groups$names, groups$names, rownames(x))
                output[[nm]] <- current
            }
        }

    } else {
        output <- do.call(score_markers_summary, c(list(x), args))
        dimnames(output$mean) <- dimnames(output$detected) <- list(rownames(x), groups$names)

        for (nm in c("cohens.d", "delta.mean", "delta.detected", "auc")) {
            current <- output[[nm]]
            if (length(current)) {
                names(current) <- groups$names
                for (i in seq_along(current)) {
                    df <- data.frame(current[[i]])
                    rownames(df) <- rownames(x)
                    current[[i]] <- df
                }
                output[[nm]] <- current
            }
        }
    }

    if (!compute.auc) {
        output <- output[setdiff(names(output), "auc")]
    }
    output
}
