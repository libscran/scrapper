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
#' @param multilevel.resolution Numeric scalar specifying the resolution when \code{method="multilevel"}.
#' Lower values favor fewer, larger communities; higher values favor more, smaller communities.
#' @param leiden.resolution Numeric scalar specifying the resolution when \code{method="leiden"}.
#' Lower values favor fewer, larger communities; higher values favor more, smaller communities.
#' @param leiden.objective String specifying the objective function when \code{method="leiden"}.
#' \code{"modularity"} uses the generalized modularity, \code{"cpm"} uses the Constant Potts Model, and \code{"er"} uses the Erd\"os-R\'enyi G(n, p) model.
#' The CPM typically yields more fine-grained clusters than the modularity at the same \code{leiden.resolution}.
#' @param walktrap.steps Integer scalar specifying the number of steps to use when \code{method="walktrap"}.
#' This determines the ability of the Walktrap algorithm to distinguish highly interconnected communities from the rest of the graph.
#' @param seed Integer scalar specifying the random seed to use for \code{method="multilevel"} or \code{"leiden"}.
#'
#' @return A list containing \code{membership}, a factor containing the cluster assignment for each cell.
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
#' The various \code{cluster_*} functions in \url{https://libscran.github.io/scran_graph_cluster/}. 
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
    leiden.objective=c("modularity", "cpm", "er"),
    walktrap.steps=4,
    seed=42)
{
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
        out <- cluster_multilevel(x, resolution=multilevel.resolution, seed=seed)
        out$levels <- lapply(out$levels, function(x) factor(x + 1L))
    } else if (method == "leiden") {
        out <- cluster_leiden(x, objective=match.arg(leiden.objective), resolution=leiden.resolution, seed=seed)
    } else if (method == "walktrap") {
        out <- cluster_walktrap(x, steps=walktrap.steps)
        out$merges <- out$merges + 1L
    }

    out$membership <- factor(out$membership + 1L)
    out
}
