#' K-means clustering
#'
#' Perform k-means clustering with a variety of different initialization and refinement algorithms.
#'
#' @param x Matrix-like object where rows are dimensions and columns are cells.
#' This is typically a dense double-precision matrix containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#' However, any matrix representation supported by \code{\link[beachmat]{initializeCpp}} can also be used.
#' @param k Integer scalar specifying the number of clusters.
#' @param init.method String specifying the initialization method for the centers:
#' \itemize{
#' \item \code{"var-part"} uses variance partitioning as described by Su and Dy (2007).
#' The dataset is repeatedly split along the dimension of greatest variance until \code{k} partitions are formed, the centroids of which form the initial clusters.
#' This approach is slower than the others but fully deterministic.
#' \item \code{"kmeans++"} uses the weighted sampling method described by Arthur and Vassilvitskii (2007).
#' \code{k} points are sampled with probability based on the smallest distance to any previously sampled point.
#' This improves the likelihood of choosing initial centroids that are far apart from each other.
#' \item \code{"random"} initialization involves choosing \code{k} random points as the initial centers.
#' This is the simplest and fastest method but may not yield good starting points. 
#' }
#' @param refine.method String specifying the refinement method.
#' \itemize{
#' \item \code{"lloyd"} uses Lloyd's algorithm, which performs a batch update in each iteration.
#' This is simple and amenable to parallelization but may not converge.
#' \item \code{"hartigan-wong"} uses the Hartigan-Wong algorithm, which transfers points between clusters to optimize the drop in the within-cluster sum of squares.
#' This is slower but has a greater chance of convergence.
#' }
#' @param var.part.optimize.partition Boolean indicating whether each partition boundary should be optimized in \code{init.method = "var.part"}.
#' This reduces the sum of squares in the child partitions, which improves the quality of the partition at the cost of some extra compute time.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param var.part.size.adjustment Number specifying the cluster size adjustment when selecting the next cluster to partition in \code{init.method = "var.part"}.
#' This should be non-negative and no greater than 1.
#' Setting it to 0 will select the cluster with the highest variance, ignoring the cluster size.
#' Setting it to 1 will select the cluster with the highest sum of squares, favoring larger clusters.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param lloyd.iterations Integer specifying the maximum number of iterations for \code{refine.method = "lloyd"}.
#' Larger values increase the chance of convergence at the cost of increasing compute time.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param hartigan.wong.iterations Integer specifying the maximum number of iterations for \code{refine.method = "hartigan-wong"}.
#' Larger values increase the chance of convergence at the cost of increasing compute time.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param hartigan.wong.quick.transfer.iterations Integer specifying the maximum number of quick transfer iterations for \code{refine.method = "hartigan-wong"}.
#' Larger values increase the chance of convergence at the cost of increasing compute time.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param hartigan.wong.quit.quick.transfer.failure Boolean indicating whether to quit kon convergence failure during quick transfer iterations in \code{refine.method = "hartigan-wong"}.
#' Setting this to \code{FALSE} gives the algorithm another chance to converge by attempting another optimal transfer iteration, at the cost of more compute time.
#' If \code{TRUE}, the function follows the same behavior as R's \code{\link{kmeans}}.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param seed Integer specifying the seed for random number generation in some of the initialization methods.
#' @param random.seed Integer specifying the seed for random number generation when \code{init.method = "random"}.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param random.seed Integer specifying the seed for random number generation when \code{init.method = "kmeans++"}.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' @param warn Boolean specifying whether a warning should be emitted if the k-means algorithm failed to converge.
#' @param num.threads Integer specifying the number of threads to use.
#' If \code{NULL}, the default value in \code{\link{clusterKmeansDefaults}} is used.
#' 
#' @return List containing:
#' \itemize{
#' \item \code{clusters}, a factor containing the cluster assignment for each cell.
#' The number of levels is no greater than \code{k}, where each level is an integer that refer to a column of \code{centers}.
#' \item \code{centers}, a numeric matrix with the coordinates of the cluster centroids (dimensions in rows, centers in columns).
#' The number of columns is no greater than \code{k}.
#' Empty clusters are automatically removed.
#' \item \code{iterations}, an integer scalar specifying the number of refinement iterations that were performed.
#' \item \code{status}, an integer scalar specifying the completion status of the algorithm.
#' A value of zero indicates success while the meaning of any non-zero value depends on the choice of \code{refine.method}:
#' \itemize{
#' \item For Lloyd, a value of 2 indicates convergence failure.
#' \item For Hartigan-Wong, a value of 2 indicates convergence failure in the optimal transfer iterations.
#' A value of 4 indicates convergence failure in the quick transfer iterations when \code{hartigan.wong.quit.quick.transfer.failure = TRUE}.
#' }
#' }
#'
#' @seealso
#' \url{https://libscran.github.io/kmeans/}, which describes the various initialization and refinement algorithms in more detail.
#'
#' \code{\link{clusterKmeans.se}}, for k-means clustering on a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @references
#' Hartigan JA. and Wong MA (1979).
#' Algorithm AS 136: A K-means clustering algorithm.
#' \emph{Applied Statistics} 28, 100-108.
#'
#' Arthur D and Vassilvitskii S (2007).
#' k-means++: the advantages of careful seeding.
#' \emph{Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms} 1027-1035.
#'
#' Su T and Dy JG (2007).
#' In Search of Deterministic Methods for Initializing K-Means and Gaussian Mixture Clustering.
#' \emph{Intelligent Data Analysis} 11, 319-338.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' clustering <- clusterKmeans(x, k=3)
#' table(clustering$clusters, iris[,"Species"])
#'
#' @export
#' @importFrom beachmat initializeCpp
clusterKmeans <- function(
    x,
    k,
    init.method = c("var-part", "kmeans++", "random"),
    refine.method = c("hartigan-wong", "lloyd"),
    warn = TRUE,
    seed = NULL,
    random.seed = seed,
    kmeanspp.seed = seed,
    var.part.optimize.partition = NULL,
    var.part.size.adjustment = NULL,
    lloyd.iterations = NULL,
    hartigan.wong.iterations = NULL,
    hartigan.wong.quick.transfer.iterations = NULL,
    hartigan.wong.quit.quick.transfer.failure = NULL,
    num.threads = NULL
) {
    .checkSEX(x, "clusterKmeans.se")

    tatami <- FALSE
    x0 <- x
    if (!is.matrix(x) || !is.numeric(x)) {
        x0 <- initializeCpp(x)
        tatami <- TRUE
    }

    output <- cluster_kmeans(
        x0,
        k,
        tatami=tatami,
        init_method=match.arg(init.method),
        refine_method=match.arg(refine.method),
        random_seed=random.seed, 
        kmeanspp_seed=kmeanspp.seed, 
        var_part_optimize_partition=var.part.optimize.partition,
        var_part_size_adjustment=var.part.size.adjustment, 
        lloyd_iterations=lloyd.iterations,
        hartigan_wong_iterations=hartigan.wong.iterations,
        hartigan_wong_quick_transfer_iterations=hartigan.wong.quick.transfer.iterations,
        hartigan_wong_quit_quick_transfer_failure=hartigan.wong.quit.quick.transfer.failure,
        nthreads=num.threads
    )

    if (output$status != 0 && warn) {
        warning("convergence failure for k-means")
    }

    output$clusters <- factor(output$clusters + 1L)
    rownames(output$centers) <- rownames(x)

    output 
}

#' Default parameters for \code{\link{clusterKmeans}}
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' clusterKmeansDefaults()
#' @export
clusterKmeansDefaults <- function() {
    def <- cluster_kmeans_defaults()
    def[["init.method"]] <- "var-part"
    def[["refine.method"]] <- "hartigan-wong"
    def[["warn"]] <- TRUE
    def
}
