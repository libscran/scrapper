#' Aggregate expression across genes
#'
#' Aggregate expression values across genes, potentially with weights.
#' This is typically used to summarize expression values for gene sets into a single per-cell score.
#'
#' @param x A matrix-like object where rows correspond to genes or genomic features and columns correspond to cells.
#' Values are usually normalized expression values, possibly log-transformed depending on the application.
#' @param sets List of vectors where each entry corresponds to a gene set.
#' Each entry may be an integer vector of row indices, a logical vector of length equal to the number of rows, or a character vector of row names.
#' For integer and character vectors, duplicate elements are ignored.
#'
#' Alternatively, each entry may be a list of two vectors.
#' The first vector should be either integer (row indices) or character (row names), specifying the genes in the set. 
#' The second vector should be numeric and of the same length as the first vector, specifying the weight associated with each gene.
#' If duplicate genes are present, only the first occurrence is used.
#' If the first vector contains gene names not present in \code{x}, those genes are ignored.
#' @param average Logical scalar indicating whether to compute the average rather than the sum.
#' @param convert Logical scalar indicating whether to convert gene identities to non-duplicate row indices in each entry of \code{sets}.
#' Can be set to \code{FALSE} for greater efficiency if the \code{sets} already contains non-duplicated integer vectors. 
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list of length equal to that of \code{sets}.
#' Each entry is a numeric vector of length equal to the number of columns in \code{x}, 
#' containing the (weighted) sum/mean of expression values for the corresponding set across all cells.
#'
#' @author Aaron Lun
#' @seealso
#' The \code{aggregate_across_genes} function in \url{https://libscran.github.io/scran_aggregate/}. 
#'
#' \code{\link{aggregateAcrossGenes.se}}, to perform aggregation on a \link[SummarizedExperiment]{SummarizedExperiment}.
#'
#' \code{\link{aggregateAcrossCells}}, to aggregate expression values across groups of cells.
#'
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' x <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Unweighted aggregation:
#' sets <- list(
#'    foo = sample(nrow(x), 20),
#'    bar = sample(nrow(x), 10)
#' )
#' agg <- aggregateAcrossGenes(x, sets)
#' str(agg)
#'
#' # Weighted aggregation:
#' sets <- list(
#'    foo = list(sample(nrow(x), 20), runif(20)),
#'    bar = list(sample(nrow(x), 10), runif(10))
#' )
#' agg2 <- aggregateAcrossGenes(x, sets, average = TRUE)
#' str(agg2)
#' 
#' @export
#' @importFrom beachmat initializeCpp
aggregateAcrossGenes <- function(x, sets, average = FALSE, convert = TRUE, num.threads = 1) {
    if (convert) {
        nr <- nrow(x)
        rn <- rownames(x)
        for (i in seq_along(sets)) {
            curset <- sets[[i]]
            if (!is.list(curset)) {
                sets[[i]] <- .sanitizeGeneSet(curset, n=nr, names=rn, arg=sprintf("sets[[%i]]", i))
            } else {
                if (length(curset) != 2L) {
                    stop("expected 'sets[[", i, "]]' to contain two vectors")
                }
                sets[[i]] <- .sanitizeGeneSetWithWeights(curset[[1]], curset[[2]], n=nr, names=rn, index=i)
            }
        }
    }

    ptr <- initializeCpp(x, .check.na=FALSE)
    output <- aggregate_across_genes(ptr, sets, average, num.threads)
    names(output) <- names(sets)
    output
}

.sanitizeGeneSetWithWeights <- function(genes, weights, n, names, index) {
    if (length(genes) != length(weights)) {
        stop("expected all vectors of 'sets[[", index, "]]' to have the same length")
    }

    if (is.numeric(genes)) {
        genes <- as.integer(genes)
        if (anyDuplicated(genes)) {
            keep <- !duplicated(genes)
            genes <- genes[keep]
            weights <- weights[keep]
        }
        if (anyNA(genes) || min(genes) < 1 || max(genes) > n) {
            stop("'sets[[", index, "]]' contains out-of-range indices")
        }

    } else {
        if (anyDuplicated(genes)) {
            keep <- !duplicated(genes)
            genes <- genes[keep]
            weights <- weights[keep]
        }
        genes <- match(genes, names)
        if (anyNA(genes)) {
            keep <- !is.na(genes)
            genes <- genes[keep]
            weights <- weights[keep]
        }
    }

    if (is.unsorted(genes)) {
        o <- order(genes)
        genes <- genes[o]
        weights <- weights[o]
    }
    list(genes, weights)
}
