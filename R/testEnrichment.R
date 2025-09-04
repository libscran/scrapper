#' Test for gene set enrichment 
#'
#' Perform a hypergeometric test for enrichment of gene sets in a list of interesting genes (e.g., markers).
#'
#' @param x Vector of identifiers for some interesting genes, e.g., symbols or Ensembl IDs.
#' This is usually derived from a selection of top markers, e.g., from \code{\link{scoreMarkers}}.
#' @param sets List of vectors of identifiers for the pre-defined gene sets.
#' Each inner vector corresponds to a gene set and should contain the same type of identifiers as \code{x}.
#' @param universe Vector of identifiers for the universe of genes in the dataset.
#' It is expected that \code{x} and each vector in \code{sets} is a subset of \code{universe}.
#' Alternatively, an integer scalar specifying the size of the universe.
#' @param log Logical scalar indicating whether to report log-transformed p-values.
#' This may be desirable to avoid underflow at near-zero p-values.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return Numeric vector of (log-transformed) p-values to test for significant enrichment of \code{x} in each entry of \code{sets}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{phyper}} and \url{https://libscran.github.io/phyper/},
#' which is the basis for the underlying calculation.
#'
#' @examples
#' testEnrichment(
#'     x=LETTERS[1:5], 
#'     sets=list(
#'         first=LETTERS[1:10],
#'         second=LETTERS[1:5 * 2],
#'         third=LETTERS[10:20]
#'     ),
#'     universe=LETTERS
#' )
#' 
#' @export
testEnrichment <- function(x, sets, universe=NULL, log=FALSE, num.threads=1) {
    num.sets <- length(sets)
    all.genes <- unlist(sets)
    set.sizes <- lengths(sets)
    set.ids <- rep(seq_along(sets), set.sizes)

    if (is.null(universe)) {
        universe <- length(unique(c(x, all.genes)))
    } else if (!is.numeric(universe) || length(universe) != 1) {
        x <- intersect(x, universe)
        keep <- all.genes %in% universe
        if (!all(keep)) {
            all.genes <- all.genes[keep]
            set.ids <- set.ids[keep]
            set.sizes <- tabulate(set.ids, nbins=num.sets)
        }
        universe <- length(universe)
    }

    overlap <- tabulate(set.ids[all.genes %in% x], nbins=num.sets)
    out <- test_enrichment(overlap, length(x), set.sizes, universe, log=log, num_threads=num.threads)
    names(out) <- names(sets)
    out
}
