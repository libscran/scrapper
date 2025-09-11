#' Summarize pairwise effect sizes for each group
#'
#' For each group, summarize the effect sizes for all pairwise comparisons to other groups.
#' This yields a set of summary statistics that can be used to rank marker genes for each group.
#'
#' @param effects A 3-dimensional numeric containing the effect sizes from each pairwise comparison between groups.
#' The extents of the first two dimensions are equal to the number of groups, while the extent of the final dimension is equal to the number of genes.
#' The entry \code{[i, j, k]} represents the effect size from the comparison of group \code{j} against group \code{i} for gene \code{k}.
#' See also the output of \code{\link{scoreMarkers}} with \code{all.pairwise=TRUE}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return List of data frames containing summary statistics for the effect sizes.
#' Each data frame corresponds to a group, each row corresponds to a gene, and each column contains a summary statistic.
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
summarizeEffects <- function(effects, num.threads=1) {
    dm <- dim(effects)
    ngroups <- dm[1]
    ngenes <- dm[3]

    out <- summarize_effects(ngenes, ngroups, effects, num_threads=num.threads)
    dmn <- dimnames(effects)
    names(out) <- dmn[[1]]

    for (i in seq_along(out)) {
        df <- data.frame(out[[i]])
        rownames(df) <- dmn[[3]]
        out[[i]] <- df
    }

    out
}
