#' Graph-based clustering of cells
#'
#' Identify clusters of cells using a variety of community detection methods from a graph where similar cells are connected.
#'
#' @param x List containing graph information or an external pointer to a graph, as returned by \code{\link{buildSnnGraph}}.
#' Alternatively, an \link[igraph]{igraph} object with edge weights.
#' @param method String specifying the algorithm to use.
#' @param multilevel.resolution Numeric scalar specifying the resolution when \code{method="multilevel"}.
#' @param leiden.resolution Numeric scalar specifying the resolution when \code{method="leiden"}.
#' @param leiden.objective String specifying the objective function when \code{method="leiden"}.
#' @param walktrap.steps Integer scalar specifying the number of steps to use when \code{method="walktrap"}.
#' @param seed Integer scalar specifying the random seed to use for \code{method="multilevel"} or \code{"leiden"}.
#'
#' @return A list containing \code{membership}, an integer vector containing the cluster assignment for each cell;
#' and \code{status}, an integer scalar indicating whether the algorithm completed successfully (0) or not (non-zero).
#' Additional fields may be present depending on the \code{method}:
#' \itemize{
#' \item For \code{method="multilevel"}, the \code{levels} list contains the clustering result at each level of the algorithm.
#' A \code{modularity} numeric vector also contains the modularity at each level, the highest of which corresponds to the reported \code{membership}.
#' \item For \code{method="leiden"}, a \code{quality} numeric scalar containg the quality of the partitioning.
#' \item For \code{method="walktrap"}, a \code{merges} matrix specifies the pair of cells or clusters that were merged at each step of the algorithm.
#' A \code{modularity} numeric scalar also contains the modularity of the final partitioning.
#' }
#'
#' @author Aaron Lun
#'
#' @seealso
#' \url{https://igraph.org/c/html/latest/igraph-Community.html}, for the underlying implementation of each clustering method.
#'
#' The various \code{cluster_*} functions in \url{https://libscran.github.io/scran_graph_cluster/}, for wrappers around the \pkg{igraph} code.
#'
#' @examples
#' data <- matrix(rnorm(10000), ncol=1000)
#' gout <- buildSnnGraph(data)
#' str(gout)
#'
#' str(clusterGraph(gout))
#' str(clusterGraph(gout, method="leiden"))
#' str(clusterGraph(gout, method="walktrap"))
#'
#' @export
#' @importFrom methods is
clusterGraph <- function(
    x,
    method=c("multilevel", "leiden", "walktrap"),
    multilevel.resolution=1, 
    leiden.resolution=1, 
    leiden.objective=c("modularity", "cpm"),
    walktrap.steps=4,
    seed=42)
{
    if (!is(x, "igraph")) {
        g <- igraph::make_undirected_graph(x$edges, n=x$vertices)
        igraph::E(g)$weight <- x$weights
        x <- g
    }

    # Some shenanigans to transparently set the seed as in the C igraph library.
    if (exists(".Random.seed")) {
        old.seed <- .Random.seed
        on.exit(set.seed(old.seed))
    }
    set.seed(seed)

    # For the time being, we still need to use the igraph R package. This will
    # change once we can switch to a vendored copy of the igraph C library,
    # which will eliminate any susceptibility to the end-user R environment.
    method <- match.arg(method)
    if (method == "multilevel") {
        out <- igraph::cluster_louvain(x, resolution=multilevel.resolution, weights=igraph::E(x)$weight)

    } else if (method == "leiden") {
        leiden.objective <- if (match.arg(leiden.objective) == "cpm") "CPM" else "modularity"
        out <- igraph::cluster_leiden(x, objective_function=leiden.objective, resolution=leiden.resolution, weights=igraph::E(x)$weight)

    } else if (method == "walktrap") {
        out <- igraph::cluster_walktrap(x, steps=walktrap.steps, weights=igraph::E(x)$weight)
    }

    out
}
