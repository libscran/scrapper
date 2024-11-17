#' Summarize pairwise effect sizes for each group
#'
#' For each group, summarize the effect sizes for all pairwise comparisons to other groups.
#' This yields a set of summary statistics that can be used to rank marker genes for each group.
#'
#' @param effects A 3-dimensional numeric containing the effect sizes from each pairwise comparison between groups.
#' The extents of the first two dimensions are equal to the number of groups, while the extent of the final dimension is equal to the number of genes.
#' The entry \code{[i, j, k]} represents Cohen's d from the comparison of group \code{j} over group \code{i} for gene \code{k}.
#' See also the output of \code{\link{scoreMarkers}} with \code{all.pairwise=TRUE}.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return List of data frames containing summary statistics for the effect sizes.
#' Each data frame corresponds to a group, each row corresponds to a gene, and each column contains a single summary.
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
#' \url{https://libscran.github.io/scran_markers/}, for more details on the statistics.
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
