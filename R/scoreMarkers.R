#' Score marker genes
#'
#' Score marker genes for each group using a variety of effect sizes from pairwise comparisons between groups.
#' This includes Cohen's d, the area under the curve (AUC), the difference in the means (delta-mean) and the difference in the proportion of detected cells (delta-detected).
#' For each group, the strongest markers are those genes with the largest effect sizes (i.e., upregulated) when compared to all other groups.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' It is typically expected to contain log-expression values, e.g., from \code{\link{normalizeCounts}}.
#' @param groups A vector specifying the group assignment for each cell in \code{x}.
#' @param block Factor specifying the block of origin (e.g., batch, sample) for each cell in \code{x}.
#' If provided, comparisons are performed within each block to ensure that block effects do not confound the estimates.
#' The weighted average of the effect sizes across all blocks is reported for each gene.
#' Alternatively \code{NULL}, if all cells are from the same block.
#' @param block.average.policy String specifying the policy to use for average statistics across blocks.
#' This can either be a (weighted) \code{"mean"} or a \code{"quantile"}.
#' Only used if \code{block} is not \code{NULL}. 
#' @param block.weight.policy String specifying the policy to use for weighting different blocks when computing the average for each statistic.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL}.
#' @param variable.block.weight Numeric vector of length 2, specifying the parameters for variable block weighting.
#' See the argument of the same name in \code{\link{computeBlockWeights}} for more detail.
#' Only used if \code{block} is not \code{NULL} and \code{block.weight.policy = "variable"}.
#' @param block.quantile Number specifying the probability of the quantile of statistics across blocks. 
#' Defaults to 0.5, i.e., the median of per-block statistics.
#' Only used if \code{block} is not \code{NULL} and \code{block.average.policy="quantile"}.
#' @param threshold Non-negative numeric scalar specifying the minimum threshold on the differences in means (i.e., the log-fold change, if \code{x} contains log-expression values). 
#' This is incorporated into the effect sizes for Cohen's d and the AUC.
#' Larger thresholds will favor genes with large differences at the expense of genes with low variance that would otherwise have comparable effect sizes.
#' @param compute.group.mean Logical scalar indicating whether to compute the group-wise mean expression for each gene.
#' @param compute.group.detected Logical scalar indicating whether to compute the group-wise proportion of detected cells for each gene.
#' @param compute.delta.mean Logical scalar indicating whether to compute the delta-means, i.e., the log-fold change when \code{x} contains log-expression values.
#' @param compute.delta.detected Logical scalar indicating whether to compute the delta-detected, i.e., differences in the proportion of cells with detected expression.
#' @param compute.cohens.d Logical scalar indicating whether to compute Cohen's d.
#' @param compute.auc Logical scalar indicating whether to compute the AUC.
#' Setting this to \code{FALSE} can improve speed and memory efficiency.
#' @param compute.summary.min Boolean specifying whether to compute the minimum as a summary statistic for each effect size.
#' Only used if \code{all.pairwise=FALSE}.
#' @param compute.summary.mean Boolean specifying whether to compute the mean as a summary statistic for each effect size.
#' Only used if \code{all.pairwise=FALSE}.
#' @param compute.summary.median Boolean specifying whether to compute the median as a summary statistic for each effect size.
#' Only used if \code{all.pairwise=FALSE}.
#' @param compute.summary.max Boolean specifying whether to compute the maximum as a summary statistic for each effect size.
#' Only used if \code{all.pairwise=FALSE}.
#' @param compute.summary.quantiles Numeric scalars containing the probabilities of quantiles to compute as summary statistics for each effect size.
#' If \code{NULL}, no quantiles are computed.
#' Only used if \code{all.pairwise=FALSE}.
#' @param compute.summary.min.rank Boolean specifying whether to compute the mininum rank as a summary statistic for each effect size.
#' If \code{NULL}, no quantiles are computed.
#' Only used if \code{all.pairwise=FALSE}.
#' @param all.pairwise Logical scalar indicating whether to report the effect sizes for every pairwise comparison between groups.
#' Alternatively, an integer scalar indicating the number of top markers to report from each pairwise comparison between groups.
#' If \code{FALSE}, only the summary statistics are reported.
#' @param min.rank.limit Integer scalar specifying the maximum value of the min-rank to report.
#' Lower values improve memory efficiency at the cost of discarding information about lower-ranked genes.
#' Only used if \code{all.pairwise=FALSE}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return
#' A named list containing:
#' \itemize{
#' \item \code{nrow}, integer specifying the number of rows in \code{x}.
#' \item \code{row.names}, character vector or \code{NULL} containing the row names of \code{x}.
#' \item \code{group.ids}, vector contaning the identities of the unique groups.
#' \item \code{mean}, a numeric matrix containing the mean expression for each group.
#' Each row is a gene and each column is a group in \code{group.ids}.
#' Omitted if \code{compute.group.mean=FALSE}.
#' \item \code{detected}, a numeric matrix containing the proportion of detected cells in each group.
#' Each row is a gene and each column is a group in \code{group.ids}.
#' Omitted if \code{compute.group.detected=FALSE}.
#' }
#' 
#' If \code{all.pairwise=FALSE}, the list also contains:
#' \itemize{
#' \item \code{cohens.d}, a list of \link[S4Vectors]{DataFrame}s where each DataFrame corresponds to a group in \code{group.ids}.
#' Each row of a DataFrame represents a gene, while each column contains a summary of Cohen's d from pairwise comparisons to all other groups.
#' This includes \code{min}, \code{mean}, \code{median}, \code{max}, \code{quantile.*} and \code{min.rank} - check out \code{?\link{summarizeEffects}} for details.
#' Omitted if \code{compute.cohens.d=FALSE}.
#' \item \code{auc}, a list like \code{cohens.d} but containing the summaries of the AUCs from each pairwise comparison.
#' Omitted if \code{compute.auc=FALSE}.
#' \item \code{delta.mean}, a list like \code{cohens.d} but containing the summaries of the delta-mean from each pairwise comparison.
#' Omitted if \code{compute.delta.mean=FALSE}.
#' \item \code{delta.detected}, a list like \code{cohens.d} but containing the summaries of the delta-detected from each pairwise comparison.
#' Omitted if \code{compute.delta.detected=FALSE}.
#' }
#'
#' If \code{all.pairwise=TRUE}, the list also contains:
#' \itemize{
#' \item \code{cohens.d}, a 3-dimensional numeric array containing the Cohen's d from each pairwise comparison between groups.
#' The extents of the first two dimensions are equal to the number of groups in \code{group.ids}, while the extent of the final dimension is equal to the number of genes.
#' The entry \code{cohens.d[i, j, k]} represents Cohen's d from the comparison of group \code{group.ids[j]} over group \code{group.ids[i]} for gene \code{k}.
#' Omitted if \code{compute.cohens.d=FALSE}.
#' \item \code{auc}, an array like \code{cohens.d} but containing the AUCs from each pairwise comparison.
#' Omitted if \code{compute.auc=FALSE}.
#' \item \code{delta.mean}, an array like \code{cohens.d} but containing the delta-mean from each pairwise comparison.
#' Omitted if \code{compute.delta.mean=FALSE}.
#' \item \code{delta.detected}, an array like \code{cohens.d} but containing the delta-detected from each pairwise comparison.
#' Omitted if \code{compute.delta.detected=FALSE}.
#' }
#'
#' If \code{all.pairwise} is an integer, the list also contains:
#' \itemize{
#' \item \code{cohens.d}, a list of list of \link[S4Vectors]{DataFrame}s containing the top genes with the largest Cohen's d for each pairwise comparison.
#' Specifically, \code{cohens.d[[i]][[j]]} is a DataFrame that contains the top \code{all.pairwise} genes from the comparison of group \code{group.ids[i]} over group \code{group.ids[j]}. 
#' Each DataFrame contains an \code{index} column, the row index of the gene; and an \code{effect} column, the Cohen's d for that gene.
#' Omitted if \code{compute.cohens.d=FALSE}.
#' \item \code{auc}, a list of list of DataFrames like \code{cohens.d} but containing the AUCs from each pairwise comparison.
#' Omitted if \code{compute.auc=FALSE}.
#' \item \code{delta.mean}, a list of list of DataFrames like \code{cohens.d} but containing the delta-mean from each pairwise comparison.
#' Omitted if \code{compute.delta.mean=FALSE}.
#' \item \code{delta.detected}, a list of list of DataFrames like \code{cohens.d} but containing the delta-detected from each pairwise comparison.
#' Omitted if \code{compute.delta.detected=FALSE}.
#' }
#'
#' All returned lists will also contain:
#'
#' @section Choice of effect size:
#' The delta-mean is the difference in the mean expression between groups.
#' This is fairly straightforward to interpret - a positive delta-mean corresponds to increased expression in the first group compared to the second. 
#' The delta-mean can also be treated as the log-fold change if the input matrix contains log-transformed normalized expression values.
#' 
#' The delta-detected is the difference in the proportion of cells with detected expression between groups.
#' This lies between 1 and -1, with the extremes occurring when a gene is silent in one group and detected in all cells of the other group.
#' For this interpretation, we assume that the input matrix contains non-negative expression values, where a value of zero corresponds to lack of detectable expression.
#' 
#' Cohen's d is the standardized difference between two groups.
#' This is defined as the difference in the mean for each group scaled by the average standard deviation across the two groups.
#' (Technically, we should use the pooled variance; however, this introduces some unintuitive asymmetry depending on the variance of the larger group, so we take a simple average instead.)
#' A positive value indicates that the gene has increased expression in the first group compared to the second.
#' Cohen's d is analogous to the t-statistic in a two-sample t-test and avoids spuriously large effect sizes from comparisons between highly variable groups.
#' We can also interpret Cohen's d as the number of standard deviations between the two group means.
#' 
#' The area under the curve (AUC) is the probability that a randomly chosen observation in one group is greater than a randomly chosen observation in the other group. 
#' Values greater than 0.5 indicate that a gene is upregulated in the first group.
#' The AUC is closely related to the U-statistic used in the Wilcoxon rank sum test. 
#' The key difference between the AUC and Cohen's d is that the former is less sensitive to the variance within each group, e.g.,
#' if two distributions exhibit no overlap, the AUC is the same regardless of the variance of each distribution. 
#' This may or may not be desirable as it improves robustness to outliers but reduces the information available to obtain a fine-grained ranking. 
#' 
#' @section With a minimum change threshold:
#' Setting a minimum change threshold (i.e., \code{threshold}) prioritizes genes with large shifts in expression instead of those with low variances.
#' Currently, only positive thresholds are supported, which focuses on genes that are upregulated in the first group compared to the second.
#' The effect size definitions are generalized when testing against a non-zero threshold:
#' \itemize{
#' \item Cohen's d is redefined as the standardized difference between the difference in means and the specified threshold,
#' analogous to the TREAT method from the \pkg{limma} package.
#' Large positive values are only obtained when the observed difference in means is significantly greater than the threshold.
#' For example, if we had a threshold of 2 and we obtained a Cohen's d of 3, this means that the observed difference in means was 3 standard deviations greater than 2.
#' Note that a negative Cohen's d cannot be intepreted as downregulation, as the difference in means may still be positive but less than the threshold.
#' \item The AUC is generalized to the probability of obtaining a random observation in one group that is greater than a random observation plus the threshold in the other group.
#' For example, if we had a threshold of 2 and we obtained an AUC of 0.8, this means that, 80% of the time,
#' the random observation from the first group would be greater than a random observation from the second group by 2 or more.
#' Again, AUCs below 0.5 cannot be interpreted as downregulation, as it may be caused by a positive shift that is less than the threshold.
#' }
#' 
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#' normed <- normalizeCounts(x, size.factors=centerSizeFactors(colSums(x)))
#'
#' # Compute marker summaries for each group:
#' g <- sample(letters[1:4], ncol(x), replace=TRUE)
#' scores <- scoreMarkers(normed, g)
#' names(scores)
#' head(scores$mean)
#' head(scores$cohens.d[["a"]])
#'
#' # Report marker statistics for a single group:
#' reportGroupMarkerStatistics(scores, "b")
#'
#' @seealso
#' The \code{score_markers_summary}, \code{score_markers_pairwise} and \code{score_markers_best} functions in \url{https://libscran.github.io/scran_markers/}.
#' See their blocked equivalents (e.g., \code{score_markers_summary_blocked}) when \code{block} is specified.
#'
#' \code{\link{summarizeEffects}}, to summarize the pairwise effects returned when \code{all.pairwise=TRUE}.
#'
#' \code{\link{reportGroupMarkerStatistics}}, to consolidate the statistics for a single group into its own data frame.
#'
#' \code{\link{scoreMarkers.se}}, to score markers from a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' @export
#' @importFrom beachmat initializeCpp
#' @importFrom S4Vectors DataFrame
scoreMarkers <- function(
    x, 
    groups, 
    block=NULL, 
    block.average.policy=c("mean", "quantile"),
    block.weight.policy=c("variable", "equal", "none"),
    variable.block.weight=c(0, 1000),
    block.quantile=0.5,
    compute.group.mean=TRUE,
    compute.group.detected=TRUE,
    compute.delta.mean=TRUE,
    compute.delta.detected=TRUE,
    compute.cohens.d=TRUE,
    compute.auc=TRUE,
    compute.summary.min=TRUE,
    compute.summary.mean=TRUE,
    compute.summary.median=TRUE,
    compute.summary.max=TRUE,
    compute.summary.quantiles=NULL,
    compute.summary.min.rank=TRUE,
    threshold=0, 
    all.pairwise=FALSE, 
    min.rank.limit=500,
    num.threads=1
) {
    .checkSEX(x, "scoreMarkers.se")

    rn <- rownames(x)
    ngenes <- nrow(x)
    x <- initializeCpp(x, .check.na=FALSE)
    groups <- .transformFactor(groups)
    block <- .transformFactor(block)

    args <- list(
        groups=groups$index,
        num_groups=length(groups$names),
        block=block$index,
        block_average_policy=match.arg(block.average.policy),
        block_weight_policy=match.arg(block.weight.policy),
        variable_block_weight=variable.block.weight,
        block_quantile=block.quantile,
        threshold=threshold,
        compute_group_mean=compute.group.mean,
        compute_group_detected=compute.group.detected,
        compute_cohens_d=compute.cohens.d,
        compute_delta_mean=compute.delta.mean,
        compute_delta_detected=compute.delta.detected,
        compute_auc=compute.auc,
        num_threads=num.threads
    )

    if (isTRUE(all.pairwise)) {
        output <- do.call(score_markers_pairwise, c(list(x), args))
        for (nm in setdiff(names(output), c("mean", "detected"))) {
            current <- output[[nm]]
            dimnames(current) <- list(groups$names, groups$names, rn)
            output[[nm]] <- current
        }

    } else if (isFALSE(all.pairwise)) {
        if (!is.null(compute.summary.quantiles)) {
            compute.summary.quantiles <- sort(unique(compute.summary.quantiles))
        }

        output <- do.call(
            score_markers_summary,
            c(
                list(
                    x,
                    min_rank_limit=min.rank.limit,
                    compute_summary_mean=compute.summary.mean,
                    compute_summary_min=compute.summary.min,
                    compute_summary_median=compute.summary.median,
                    compute_summary_max=compute.summary.max,
                    compute_summary_quantiles=compute.summary.quantiles,
                    compute_summary_min_rank=compute.summary.min.rank
                ),
                args
            )
        )

        for (nm in setdiff(names(output), c("mean", "detected"))) {
            current <- output[[nm]]
            names(current) <- groups$names
            output[[nm]] <- format_summary_output(current, ngenes, rn, compute.summary.quantiles)
        }

    } else {
        output <- do.call(score_markers_best, c(list(x, top=all.pairwise), args))
        for (nm in setdiff(names(output), c("mean", "detected"))) {
            current <- output[[nm]]
            names(current) <- groups$names
            for (i in seq_along(current)) {
                names(current[[i]]) <- groups$names
            }
            output[[nm]] <- current
        }
    }

    output$nrow <- ngenes
    output$row.names <- rn
    output$group.ids <- groups$names

    if (compute.group.mean) {
        dimnames(output$mean) <- list(rn, groups$names)
    }
    if (compute.group.detected) {
        dimnames(output$detected) <- list(rn, groups$names)
    }

    output
}

#' Report marker statistics for a single group
#'
#' Combine all marker statistics for a single group into a data frame for easy inspection.
#' Users can pick one of the columns for sorting potential marker genes. 
#'
#' @param results Named list of marker statistics, typically generated by \code{\link{scoreMarkers}} with \code{all.pairwise=FALSE}.
#' @param group String or integer scalar specifying the group of interest.
#' @param effect.sizes Character vector specifying the effect sizes of interest.
#' If \code{NULL}, all effect sizes are reported in the returned data frame.
#' @param summaries Character vector specifying the summary statistics of interest.
#' If \code{NULL}, all summaries are reported in the returned data frame.
#' @param include.mean Logical scalar indicating whether the mean expression should be reported in the returned data frame.
#' @param include.detected Logical scalar indicating whether the proportion of detected cells should be reported in the returned data frame.
#'
#' @return Data frame where each row corresponds to a gene.
#' Each column contains the requested statistics for \code{group}.
#' Effect size summary columns are named as \code{<EFFECT>.<SUMMARY>}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{scoreMarkers}}, to generate \code{results}.
#'
#' \code{\link{summarizeEffects}}, for the trade-offs between effect size summaries.
#'
#' @export
reportGroupMarkerStatistics <- function(
    results,
    group,
    effect.sizes = NULL,
    summaries = NULL,
    include.mean = TRUE,
    include.detected = TRUE)
{
    if (is.null(effect.sizes)) {
        effect.sizes <- c("cohens.d", "auc", "delta.mean", "delta.detected")
    }
    if (is.null(summaries)) {
        summaries <- c("min", "mean", "median", "max", "min.rank")
    }

    current <- data.frame(row.names=rownames(results$mean), matrix(0L, nrow(results$mean), 0L))

    if (include.mean) {
        current$mean <- results$mean[,group]
    }
    if (include.detected) {
        current$detected <- results$detected[,group]
    }

    for (eff in effect.sizes) {
        eff.all <- results[[eff]]
        if (is.null(eff.all)) {
            next
        }
        eff.df <- eff.all[[group]]
        for (summ in summaries) {
            current[[paste0(eff, ".", summ)]] <- eff.df[[summ]]
        }
    }

    current
}
