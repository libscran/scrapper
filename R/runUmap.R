#' Uniform manifold approxation and projection
#'
#' Compute UMAP coordinates to visualize similarities between cells.
#'
#' @inheritParams runTsne
#' @param num.dim Integer scalar specifying the number of dimensions of the output embedding.
#' @param local.connectivity Numeric scalar specifying the number of nearest neighbors that are assumed to be always connected, with maximum membership confidence.
#' Larger values increase the connectivity of the embedding and reduce the focus on local structure.
#' This may be a fractional number of neighbors, in which case interpolation is performed when computing the membership confidence.
#' @param bandwidth Numeric scalar specifying the effective bandwidth of the kernel when converting the distance to a neighbor into a fuzzy set membership confidence.
#' Larger values reduce the decay in confidence with respect to distance, increasing connectivity and favoring global structure. 
#' @param mix.ratio Numeric scalar between 0 and 1 specifying the mixing ratio when combining fuzzy sets.
#' A mixing ratio of 1 will take the union of confidences, a ratio of 0 will take the intersection, and intermediate values will interpolate between them.
#' Larger values favor connectivity and more global structure.
#' @param spread Numeric scalar specifying the scale of the coordinates of the final low-dimensional embedding.
#' Ignored if \code{a} and \code{b} are provided.
#' @param min.dist Numeric scalar specifying the minimum distance between observations in the final low-dimensional embedding.
#' Smaller values will increase local clustering while larger values favor a more even distribution of observations throughout the low-dimensional space.
#' This is interpreted relative to \code{spread}.
#' Ignored if \code{a} and \code{b} are provided.
#' @param a Numeric scalar specifying the \eqn{a} parameter for the fuzzy set membership strength calculations.
#' Larger values yield a sharper decay in membership strength with increasing distance between observations.
#' If this or \code{b} are \code{NULL}, a suitable value for this parameter is automatically determined from \code{spread} and \code{min.dist}.
#' @param b Numeric scalar specifying the \eqn{b} parameter for the fuzzy set membership strength calculations.
#' Larger values yield an earlier decay in membership strength with increasing distance between observations.
#' If this or \code{a} are \code{NULL}, a suitable value for this parameter is automatically determined from \code{spread} and \code{min.dist}.
#' @param repulsion.strength Numeric scalar specifying the modifier for the repulsive force.
#' Larger values increase repulsion and favor local structure.
#' @param initialize.method String specifying how to initialize the embedding.
#' This should be one of:
#' \itemize{
#' \item \code{spectral}: spectral decomposition of the normalized graph Laplacian.
#' Specifically, the initial coordinates are defined from the eigenvectors corresponding to the smallest non-zero eigenvalues.
#' This fails in the presence of multiple graph components or if the approximate SVD fails to converge.
#' \item \code{random}: fills the embedding with random draws from a normal distribution.
#' \item \code{none}: uses initial values from \code{initial.coordinates}.
#' }
#' @param initialize.random.on.spectral.fail Logical scalar indicating whether to fall back to random sampling (i.e., same as \code{random})
#' if spectral initialization fails due to the presence of multiple components in the graph.
#' If \code{FALSE}, the values in \code{initial.coordinates} will be used instead, i.e., same as \code{none}.
#' Only relevant if \code{initialize.method = "spectral"} and spectral initialization fails.
#' @param initialize.spectral.scale Numeric scalar specifying the maximum absolute magnitude of the coordinates after spectral initialization.
#' All initial coordinates are scaled such that the maximum of the absolute values is equal to \code{initialize.spectral.scale}.
#' This ensures that outlier observations will not have large absolute distances that may interfere with optimization.
#' Only relevant if \code{initialize.method = "spectral"} and spectral initialization does not fail.
#' @param initialize.spectral.jitter Logical scalar indicating whether to jitter coordinates after spectral initialization to separate duplicate observations (e.g., to avoid overplotting).
#' This is done using normally-distributed noise of mean zero and standard deviation of \code{initialize.spectral.jitter.sd}.
#' Only relevant if \code{initialize.method = "spectral"} and spectral initialization does not fail.
#' @param initialize.spectral.jitter.sd Numeric scalar specifying the standard deviation of the jitter to apply after spectral initialization.
#' Only relevant if \code{initialize.method = "spectral"} and spectral initialization does not fail and \code{initialize.spectral.jitter = TRUE}.
#' @param initialize.random.scale Numeric scalar specifying the scale of the randomly generated initial coordinates.
#' Coordinates are sampled from a uniform distribution from \eqn{[-x, x)} where \eqn{x} is \code{initialize.random.scale}.
#' Only relevant if \code{initialize.method = "random"},
#' or \code{initialize.method = "spectral"} and spectral initialization fails and \code{initialize.random.on.spectral.fail = TRUE}.
#' @param initialize.seed Numeric scalar specifying the seed for the random number generation during initialization.
#' Only relevant if \code{initialize.method = "random"},
#' or \code{initialize.method = "spectral"} and \code{initialize.spectral.jitter = TRUE};
#' or \code{initialize.method = "spectral"} and spectral initialization fails and \code{initialize.random.on.spectral.fail = TRUE}.
#' @param initial.coordinates Numeric matrix of initial coordinates, with number of rows equal to the number of observations and number of columns equal to \code{num.dim}.
#' Only relevant if \code{initialize.method = "none"};
#' or \code{initialize.method = "spectral"} and spectral initialization fails and \code{initialize.random.on.spectral.fail = FALSE}.
#' @param num.epochs Integer scalar specifying the number of epochs for the gradient descent, i.e., optimization iterations. 
#' Larger values improve accuracy at the cost of increased compute time.
#' If \code{NULL}, a value is automatically chosen based on the size of the dataset:
#' \itemize{
#' \item For datasets with no more than 10000 observations, the default number of epochs is set to 500.
#' \item For larger datasets, the number of epochs is inversely proportional to the number of cells, starting from 500 and decreasing asymptotically to a lower limit of 200.
#' This choice aims to reduce computational work for very large datasets. 
#' }
#' @param learning.rate Numeric scalar specifying the initial learning rate used in the gradient descent.
#' Larger values can accelerate convergence but at the risk of skipping over suitable local optima.
#' @param negative.sample.rate Numeric scalar specifying the rate of sampling negative observations to compute repulsive forces.
#' Greater values will improve accuracy but increase compute time. 
#' @param num.neighbors Integer scalar specifying the number of neighbors to use to define the fuzzy sets.
#' Larger values improve connectivity and favor preservation of global structure, at the cost of increased compute time.
#' If \code{x} contains pre-computed neighbor search result, the number of neighbors should be equal to \code{num.neighbors}.
#' @param optimize.seed Numeric scalar specifying the seed to use for the optimization epochs. 
#' @param num.threads Integer scalar specifying the number of threads to use.
#' @param parallel.optimization Logical scalar specifying whether to parallelize the optimization step.
#' 
#' @return 
#' A numeric matrix where rows are cells and columns are the two dimensions of the embedding.
#' 
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://libscran.github.io/umappp/}, for details on the underlying implementation.
#'
#' \code{\link{runUmap.se}}, to run UMAP on a \link[SingleCellExperiment]{SingleCellExperiment}.
#'
#' @references
#' McInnes L, Healy J, Melville J (2020).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.
#' \emph{arXiv}, \url{https://arxiv.org/abs/1802.03426}
#'
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' embedding <- runUmap(x)
#' plot(embedding[,1], embedding[,2], col=iris[,5])
#' 
#' @export
#' @importFrom BiocNeighbors findKNN AnnoyParam
runUmap <- function(
    x, 
    num.dim=2,
    local.connectivity=1,
    bandwidth=1,
    mix.ratio=1,
    spread=1,
    min.dist=0.1, 
    a=NULL,
    b=NULL,
    repulsion.strength=1,
    initialize.method=c("spectral", "random", "none"),
    initial.coordinates=NULL,
    initialize.random.on.spectral.fail=TRUE,
    initialize.spectral.scale=10,
    initialize.spectral.jitter=FALSE,
    initialize.spectral.jitter.sd=0.0001,
    initialize.random.scale=10,
    initialize.seed=9876543210,
    num.epochs=NULL, 
    learning.rate=1,
    negative.sample.rate=5,
    num.neighbors=15, 
    optimize.seed=1234567890, 
    num.threads=1,
    parallel.optimization=FALSE,
    BNPARAM=AnnoyParam())
{
    if (!is.list(x)) {
        x <- findKNN(x, k=num.neighbors, transposed=TRUE, get.index="transposed", get.distance="transposed", num.threads=num.threads, BNPARAM=BNPARAM)
    } else {
        .checkNeighborIndices(x$index, num.neighbors)
    }

    if (!is.null(initial.coordinates)) {
        num.obs <- as.integer(ncol(x$index))
        stopifnot(identical(nrow(initial.coordinates), num.obs))
        stopifnot(identical(ncol(initial.coordinates), as.integer(num.dim)))
        initial.coordinates <- t(initial.coordinates)
    }

    output <- run_umap(
        nnidx=x$index, 
        nndist=x$distance, 
        num_dim=num.dim,
        local_connectivity=local.connectivity,
        bandwidth=bandwidth,
        mix_ratio=mix.ratio,
        spread=spread,
        min_dist=min.dist,
        a=a,
        b=b,
        repulsion_strength=repulsion.strength,
        initialize_method=match.arg(initialize.method),
        initial_coordinates=initial.coordinates,
        initialize_random_on_spectral_fail=initialize.random.on.spectral.fail,
        initialize_spectral_scale=initialize.spectral.scale,
        initialize_spectral_jitter=initialize.spectral.jitter,
        initialize_spectral_jitter_sd=initialize.spectral.jitter.sd,
        initialize_random_scale=initialize.random.scale,
        initialize_seed=initialize.seed,
        num_epochs=num.epochs,
        learning_rate=learning.rate,
        negative_sample_rate=negative.sample.rate,
        optimize_seed=optimize.seed,
        num_threads=num.threads,
        parallel_optimization=parallel.optimization
    )

    t(output)
}
