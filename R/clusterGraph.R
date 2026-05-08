#' Graph-based clustering of cells
#'
#' Identify clusters by applying community detection algorithms to a graph.
#' This assumes that that the nodes on the graph represent cells and weighted edges are formed between related cells.
#'
#' @param x List containing graph information or an external pointer to a graph, as returned by \code{\link{buildSnnGraph}}.
#' Alternatively, an \link[igraph]{igraph} object with edge weights.
#' @param method String specifying the algorithm to use.
#' \itemize{
#' \item \code{"multilevel"} uses multi-level modularity optimization, also known as the Louvain algorithm,
#' see \url{https://igraph.org/c/doc/igraph-Community.html#igraph_community_multilevel} for details.
#' \item \code{"walktrap"} uses the Walktrap community finding algorithm,
#' see \url{https://igraph.org/c/doc/igraph-Community.html#igraph_community_walktrap} for details.
#' \item \code{"leiden"} uses the Leiden algorithm,
#' see \url{https://igraph.org/c/doc/igraph-Community.html#igraph_community_leiden} for details.
#' }
#' @param multilevel.resolution Number specifying the resolution when \code{method="multilevel"}.
#' Lower values favor fewer, larger communities; higher values favor more, smaller communities.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#' @param leiden.resolution Number specifying the resolution when \code{method="leiden"}.
#' Lower values favor fewer, larger communities; higher values favor more, smaller communities.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#' @param leiden.objective String specifying the objective function when \code{method="leiden"}.
#' \code{"modularity"} uses the generalized modularity, \code{"cpm"} uses the Constant Potts Model, and \code{"er"} uses the Erd\"os-R\'enyi G(n, p) model.
#' The CPM typically yields more fine-grained clusters than the modularity at the same \code{leiden.resolution}.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#' @param walktrap.steps Integer specifying the number of steps to use when \code{method="walktrap"}.
#' This determines the ability of the Walktrap algorithm to distinguish highly interconnected communities from the rest of the graph.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#' @param seed Integer specifying the random seed for some of community detection algorithms.
#' @param multilevel.seed Integer specifying the random seed to use for \code{method="multilevel"}.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#' @param leiden.seed Integer specifying the random seed to use for \code{method="leiden"}.
#'
#' If \code{NULL}, the default value in \code{\link{clusterGraphDefaults}} is used.
#'
#' @return A list containing \code{membership}, a factor containing the cluster assignment for each cell.
#' Additional fields may be present depending on the \code{method}:
#' \itemize{
#' \item For \code{method="multilevel"}, the \code{levels} list contains the clustering result at each level of the algorithm.
#' A \code{modularity} numeric vector also contains the modularity at each level, the highest of which corresponds to the reported \code{membership}.
#' \item For \code{method="leiden"}, the \code{quality} number contains the quality of the partitioning.
#' \item For \code{method="walktrap"}, the \code{merges} matrix specifies the pair of cells or clusters that were merged at each step of the algorithm.
#' The \code{modularity} number also contains the modularity of the final partitioning.
#' }
#'
#' @author Aaron Lun
#'
#' @seealso
#' The various \code{cluster_*} functions in \url{https://libscran.github.io/scran_graph_cluster/}. 
#'
#' \code{\link{clusterGraph.se}}, which performs clustering on graph constructed from a \link[SingleCellExperiment]{SingleCellExperiment}.
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
    method = c("multilevel", "leiden", "walktrap"),
    seed = NULL,
    multilevel.resolution = NULL, 
    multilevel.seed = seed,
    leiden.resolution = NULL,
    leiden.seed = seed,
    leiden.objective = "modularity", 
    walktrap.steps = NULL
) {
    .checkSEX(x, "clusterGraph.se")
    if (is(x, "igraph")) {
        x <- list(
            vertices=igraph::vcount(x),
            edges=as.vector(t(igraph::as_edgelist(x, names=FALSE))),
            weights=igraph::E(x)$weight
        )
    }
    if (is.list(x)) {
        x <- list_to_graph(x)
    }

    method <- match.arg(method)
    if (method == "multilevel") {
        out <- cluster_multilevel(x, resolution=multilevel.resolution, seed=multilevel.seed)
        out$levels <- lapply(out$levels, function(x) factor(x + 1L))
    } else if (method == "leiden") {
        out <- cluster_leiden(x, objective=leiden.objective, resolution=leiden.resolution, seed=leiden.seed)
    } else if (method == "walktrap") {
        out <- cluster_walktrap(x, steps=walktrap.steps)
        out$merges <- out$merges + 1L
    }

    out$membership <- factor(out$membership + 1L)
    out
}

#' Default parameters for \code{\link{clusterGraph}}
#' @description Default parameters from the underlying C++ library.
#' Some of these may be overridden by defaults in the \code{\link{clusterGraph}} function signature.
#' @return Named list containing default values for various function arguments.
#' @author Aaron Lun
#' @examples
#' clusterGraphDefaults()
#' @export
clusterGraphDefaults <- function() cluster_graph_defaults()
