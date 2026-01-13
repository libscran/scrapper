#' Count cells in groups and blocks
#'
#' Tabulate the frequency of cells in each combination of groups and blocks.
#' This is typically used to examine the distribution of cells across batches for each cluster -
#' the presence of a batch-specific cluster may be indicative of a batch effect.
#'
#' @param groups Factor specifying the group to which each cell was assigned.
#' This is typically used for clusters.
#' @param block Factor specifying the block to which each cell was assigned.
#' This is typically used for batches or samples.
#' @param normalize.block Boolean indicating whether to normalize the number of cells across blocks.
#' If \code{TRUE}, frequencies are divided by the column sums.
#' @param normalize.groups Boolean indicating whether to normalize the number of cells across groups.
#' If \code{TRUE}, frequencies are divided by the row sums.
#' This is performed after normalization of the block counts if \code{normalize.block=TRUE}.
#' 
#' @return Matrix of (normalized) frequencies.
#' Each row corresponds to a group and each column corresponds to a block.
#'
#' @author Aaron Lun
#' @examples
#' groups <- sample(10, 100, replace=TRUE)
#' block <- sample(LETTERS[1:6], 100, replace=TRUE)
#'
#' countGroupsByBlock(groups, block)
#' countGroupsByBlock(groups, block, normalize.block=TRUE)
#' countGroupsByBlock(groups, block, normalize.groups=TRUE)
#' countGroupsByBlock(groups, block, normalize.block=TRUE, normalize.groups=TRUE)
#'
#' @seealso
#' \code{\link{table}}, which is used internally by this function.
#'
#' @export
countGroupsByBlock <- function(groups, block, normalize.block=FALSE, normalize.groups=FALSE) {
    tab <- table(groups, block)
    if (normalize.block) {
        tab <- t(t(tab) / colSums(tab))
    }
    if (normalize.groups) {
        tab <- tab / rowSums(tab)
    }
    tab
}
