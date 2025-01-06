#' K-means clustering
#'
#' Perform k-means clustering with a variety of different initialization and refinement algorithms.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells.
#' @param k Integer scalar specifying the number of clusters.
#' @param init.method String specifying the initialization method:
#' variance partitioning (\code{"var-part"}), kmeans++ (\code{"kmeans++"}) or random initialization (\code{"random"}).
#' @param refine.method String specifying the refinement method:
#' Lloyd's algorithm (\code{"lloyd"}) or the Hartigan-Wong algorithm (\code{"hartigan-wong"}).
#' @param var.part.optimize.partition Logical scalar indicating whether each partition boundary should be optimized to reduce the sum of squares in the child partitions.
#' Only used if \code{init.method = "var.part"}.
#' @param var.part.size.adjustment Numeric scalar between 0 and 1, specifying the adjustment to the cluster size when prioritizing the next cluster to partition.
#' Setting this to 0 will ignore the cluster size while setting this to 1 will generally favor larger clusters.
#' Only used if \code{init.method = "var.part"}.
#' @param lloyd.iterations Integer scalar specifying the maximmum number of iterations for the Lloyd algorithm.
#' @param hartigan.wong.iterations Integer scalar specifying the maximmum number of iterations for the Hartigan-Wong algorithm.
#' @param hartigan.wong.quick.transfer.iterations Integer scalar specifying the maximmum number of quick transfer iterations for the Hartigan-Wong algorithm.
#' @param hartigan.wong.quit.quick.transfer.failure Logical scalar indicating whether to quit the Hartigan-Wong algorithm upon convergence failure during quick transfer iterations.
#' @param seed Integer scalar specifying the seed to use for random or kmeans++ initialization.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @return 
#' By default, a list is returned containing:
#' \itemize{
#' \item \code{clusters}, a factor containing the cluster assignments for each cell.
#' \item \code{centers}, a numeric matrix with the coordinates of the cluster centroids (dimensions in rows, centers in columns).
#' \item \code{iterations}, an integer scalar specifying the number of refinement iterations that were performed.
#' \item \code{status}, an integer scalar specifying the convergence status.
#' Any non-zero value indicates a convergence failure though the exact meaning depends on the choice of \code{refine.method}.
#' }
#'
#' @seealso
#' \url{https://ltla.github.io/CppKmeans/}, which describes the various initialization and refinement algorithms in more detail.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterKmeans(x, k=3)
#' table(clustering$clusters, iris[,"Species"])
#'
#' @export
clusterKmeans <- function(
    x,
    k,
    init.method = c("var-part", "kmeans++", "random"),
    refine.method = c("hartigan-wong", "lloyd"),
    var.part.optimize.partition = TRUE,
    var.part.size.adjustment = 1,
    lloyd.iterations = 100,
    hartigan.wong.iterations = 10,
    hartigan.wong.quick.transfer.iterations = 50,
    hartigan.wong.quit.quick.transfer.failure = FALSE,
    seed=5489L,
    num.threads=1)
{
    output <- cluster_kmeans(
        x,
        k,
        init_method=match.arg(init.method),
        refine_method=match.arg(refine.method),
        var_part_optimize_partition=var.part.optimize.partition,
        var_part_size_adjustment=var.part.size.adjustment, 
        lloyd_iterations=lloyd.iterations,
        hartigan_wong_iterations=hartigan.wong.iterations,
        hartigan_wong_quick_transfer_iterations=hartigan.wong.quick.transfer.iterations,
        hartigan_wong_quit_quick_transfer_failure=hartigan.wong.quit.quick.transfer.failure,
        seed=seed, 
        nthreads=num.threads
    )
    output$clusters <- factor(output$clusters + 1L)
    output 
}
