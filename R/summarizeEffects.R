#' Summarize pairwise effect sizes for each group
#'
#' For each group, summarize the effect sizes for all pairwise comparisons to other groups.
#' This yields a set of summary statistics that can be used to rank marker genes for each group.
#'
#' @param effects A 3-dimensional numeric containing the effect sizes from each pairwise comparison between groups.
#' The extents of the first two dimensions are equal to the number of groups, while the extent of the final dimension is equal to the number of genes.
#' The entry \code{[i, j, k]} represents the effect size from the comparison of group \code{j} against group \code{i} for gene \code{k}.
#' See also the output of \code{\link{scoreMarkers}} with \code{all.pairwise=TRUE}.
#' @param compute.summary.min Boolean specifying whether to compute the minimum as a summary statistic.
#' @param compute.summary.mean Boolean specifying whether to compute the mean as a summary statistic. 
#' @param compute.summary.median Boolean specifying whether to compute the median as a summary statistic. 
#' @param compute.summary.max Boolean specifying whether to compute the maximum as a summary statistic. 
#' @param compute.summary.quantiles Numeric scalars containing the probabilities of quantiles to compute as summary statistics. 
#' If \code{NULL}, no quantiles are computed.
#' @param compute.summary.min.rank Boolean specifying whether to compute the mininum rank as a summary statistic. 
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return List of \link[S4Vectors]{DataFrame}s containing summary statistics for the effect sizes.
#' Each DataFrame corresponds to a group, each row corresponds to a gene, and each column contains a summary statistic.
#' If \code{compute.summary.quantiles} is provided, the \code{"quantile"} column is a nested DataFrame where each column coresponds to a probability in \code{compute.summary.quantiles}.
#'
#' @details
#' Each summary statistic can be used to prioritize different sets of marker genes for the group of interest, by ranking them in decreasing order according to said statistic:
#' \itemize{
#' \item \code{min} contains the minimum effect size across all comparisons involving the group of interest.
#' Genes with large values are upregulated in all comparisons.
#' As such, it is the most stringent summary as markers will only have large values if they are uniquely upregulated in the group of interest compared to every other group.
#' \item \code{mean} contains the mean effect size across all comparisons involving the group of interest.
#' Genes with large values are upregulated on average compared to the other groups.
#' This is a good general-purpose summary statistic.
#' \item \code{median} contains the median effect size across all comparisons involving the group of interest.
#' Genes with large values are upregulated compared to most (i.e., at least 50% of) other groups.
#' Compared to the mean, this is more robust to outlier effects but less sensitive to strong effects in a minority of comparisons.
#' \item \code{max} contains the maximum effect size across all comparisons involving the group of interest.
#' Using this to define markers will focus on genes that are upregulated in at least one comparison.
#' As such, it is the least stringent summary as markers can achieve large values if they are upregulated in the group of interest compared to any one other group.
#' \item \code{quantile[[P]]} contains the quantile P across all comparisons involving the group of interest.
#' This is a generalization of the minimum, median and maximum for arbitrary quantile probabilities.
#' For example, a large \code{quantile[["20"]]} would mean that the gene is upregulated in the group of interest compared to 80% of other groups.
#' }
#'
#' The exact definition of \dQuote{large} depends on the choice of effect size.
#' For signed effects like Cohen's d, delta-mean and delta-detected, the value must be positive to be considered \dQuote{large}.
#' For the AUC, a value greater than 0.5 is considered \dQuote{large}.
#' This interpretation is also affected by the choice of \code{threshold=} used to compute each effect size in \code{\link{scoreMarkers}},
#' e.g., a negative Cohen's d cannot be interpreted as downregulation when the threshold is positive. 
#' 
#' The \code{min.rank} is a more exotic summary statistic, containing the minimum rank for each gene across all comparisons involving the group of interest.
#' This is defined by ranking the effect sizes across genes within each comparison, and then taking the minimum of these ranks across comparisons.
#' Taking all genes with \code{min.rank <= T} will yield a set containing the top \code{T} genes from each comparison.
#' The idea is to ensure that there are at least \code{T} genes that can distinguish the group of interest from any other group.
#'
#' \code{NaN} effect sizes are allowed, e.g., if two groups do not exist in the same block for a blocked analysis in \code{\link{scoreMarkers}} with \code{block=}. 
#' This function will ignore \code{NaN} values when computing each summary.
#' If all effects are NaN for a particular group, the summary statistic will also be \code{NaN}.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' normed <- normalizeCounts(x, size.factors=centerSizeFactors(colSums(x)))
#'
#' g <- sample(letters[1:4], ncol(x), replace=TRUE)
#' effects <- scoreMarkers(normed, g, all.pairwise=TRUE)
#'
#' summarized <- summarizeEffects(effects$cohens.d)
#' str(summarized)
#'
#' @seealso
#' The \code{summarize_effects} function in \url{https://libscran.github.io/scran_markers/}.
#'
#' \code{\link{scoreMarkers}}, to compute the pairwise effects in the first place.
#'
#' @export
summarizeEffects <- function(
    effects,
    compute.summary.min=TRUE,
    compute.summary.mean=TRUE,
    compute.summary.median=TRUE,
    compute.summary.max=TRUE,
    compute.summary.quantiles=NULL,
    compute.summary.min.rank=TRUE,
    num.threads=1
) {
    dm <- dim(effects)
    ngroups <- dm[1]
    ngenes <- dm[3]

    out <- summarize_effects(
        ngenes,
        ngroups,
        effects,
        compute_summary_mean=compute.summary.mean,
        compute_summary_min=compute.summary.min,
        compute_summary_median=compute.summary.median,
        compute_summary_max=compute.summary.max,
        compute_summary_quantiles=compute.summary.quantiles,
        compute_summary_min_rank=compute.summary.min.rank,
        num_threads=num.threads
    )
    dmn <- dimnames(effects)
    names(out) <- dmn[[1]]

    format_summary_output(out, ngenes, dmn[[3]], compute.summary.quantiles)
}

#' @importFrom S4Vectors DataFrame
format_summary_output <- function(raw, ngenes, rownames, quantiles) {
    for (i in seq_along(raw)) {
        current <- raw[[i]]

        if ("quantile" %in% names(current)) {
            names(current$quantile) <- as.character(quantiles)
            to.add <- DataFrame(current$quantile, check.names=FALSE)
            current$quantile <- logical(ngenes)
        } else {
            to.add <- NULL
        }

        if (length(current)) {
            df <- DataFrame(current)
        } else {
            df <- DataFrame(matrix(0, ngenes, 0))
        }
        rownames(df) <- rownames

        if (!is.null(to.add)) {
            df$quantile <- to.add
        }

        raw[[i]] <- df
    }
    raw
}
