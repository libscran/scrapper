#' Test for gene set enrichment 
#'
#' Perform a hypergeometric test for enrichment of interesting genes (e.g., markers) in one or more pre-defined gene sets.
#'
#' @param x Vector of identifiers for the interesting genes.
#' @param sets List of vectors of identifiers for the pre-defined gene sets.
#' @param universe Vector of identifiers for the universe of genes in the dataset.
#' It is expected that \code{x} is a subset of \code{universe}.
#' Alternatively, an integer scalar specifying the size of the universe.
#' @param log Logical scalar indicating whether to report the log-transformed p-values.
#' @param num.threads Integer scalar specifying the number of threads to use.
#'
#' @return Numeric vector of (log-transformed) p-values to test for significant enrichment of \code{x} in each entry of \code{sets}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{phyper}} and \url{https://github.com/libscran/phyper/},
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
testEnrichment <- function(x, sets, universe, log=FALSE, num.threads=1) {
    all.genes <- unlist(sets)
    set.sizes <- lengths(sets)
    set.ids <- rep(seq_along(sets), set.sizes)

    overlap <- tabulate(set.ids[all.genes %in% x], nbins=length(sets))
    if (!is.numeric(universe) && length(universe) != 1) {
        set.sizes <- tabulate(set.ids[all.genes %in% universe], nbins=length(sets))
        universe <- length(universe)
    }

    out <- test_enrichment(overlap, length(x), set.sizes, universe, log=log, num_threads=num.threads)
    names(out) <- names(sets)
    out
}
