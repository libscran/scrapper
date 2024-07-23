#' Summarize pairwise effect sizes for each group
#'
#' For each group, summarize the effect sizes for all pairwise comparisons to other groups.
#' This yields a set of summary statistics that can be used to rank marker genes
#'
#' @param effects A 3-dimensional numeric containing the effect sizes from each pairwise comparison between groups.
#' The first dimension represents the first group, the second dimension represents the second group, and the final dimension represents the gene;
#' the entry \code{[i, j, k]} represents Cohen's d for \code{i} minus \code{j} of gene \code{k}.
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
#' The \code{summarize_effects} function from \url{https://libscran.github.io/scran_markers}.
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
