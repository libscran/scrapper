#' Run all neighbor-related steps
#'
#' Run all steps that require a nearest-neighbor search.
#' This includs \code{\link{runUmap}}, \code{\link{runTsne}} and \code{\link{buildSnnGraph}} with \code{\link{clusterGraph}}.
#' The idea is to build the index once, perform the neighbor search, and run each task in parallel to save time.
#'
#' @param x Numeric matrix where rows are dimensions and columns are cells,
#' typically containing a low-dimensional representation from, e.g., \code{\link{runPca}}.
#'
#' Alternatively, an index constructed by \code{\link[BiocNeighbors]{buildIndex}}.
#' @param runUmap.args Named list of further arguments to pass to \code{\link{runUmap}}.
#' This can be set to \code{NULL} to omit the UMAP.
#' @param runTsne.args Named list of further arguments to pass to \code{\link{runTsne}}.
#' This can be set to \code{NULL} to omit the t-SNE.
#' @param buildSnnGraph.args Named list of further arguments to pass to \code{\link{buildSnnGraph}}.
#' Ignored if \code{clusterGraph.args=NULL}.
#' @param clusterGraph.args Named list of further arguments to pass to \code{\link{clusterGraph}}.
#' This can be set to \code{NULL} to omit the clustering.
#' @param BNPARAM A \link[BiocNeighbors]{BiocNeighborParam} instance specifying the nearest-neighbor search algorithm to use.
#' @param collapse.search Logical scalar indicating whether to collapse the nearest-neighbor search for each step into a single search.
#' Steps that need fewer neighbors will take a subset of the neighbors from the collapsed search.
#' This is faster but may not give the same results as separate searches for some algorithms (e.g., approximate searches).
#' @param return.graph Logical scalar indicating whether to return the output of \code{\link{buildSnnGraph}}.
#' By default, only the output of \code{\link{clusterGraph}} is returned.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' At least one thread should be available for each step. 
#'
#' @return A named list containing the results of each step.
#' See each individual function for the format of the results.
#' 
#' @author Aaron Lun
#' @examples
#' x <- t(as.matrix(iris[,1:4]))
#' # (Turning down the number of threads so that R CMD check is happy.)
#' res <- runAllNeighborSteps(x, num.threads=2)
#' str(res)
#'
#' @export
#' @importFrom BiocNeighbors buildIndex AnnoyParam findKNN
#' @importFrom parallel makeCluster stopCluster parLapply
runAllNeighborSteps <- function(
    x,
    runUmap.args=list(),
    runTsne.args=list(),
    buildSnnGraph.args=list(),
    clusterGraph.args=list(),
    BNPARAM=AnnoyParam(),
    return.graph=FALSE,
    collapse.search=FALSE,
    num.threads=3) 
{
    k.choices <- list()

    if (!is.null(runUmap.args)) {
        umap.k <- runUmap.args$num.neighbors
        if (is.null(umap.k)) {
            umap.k <- formals(runUmap)$num.neighbors
        }
        k.choices$runUmap <- umap.k
    }

    if (!is.null(runTsne.args)) {
        tsne.perplexity <- runTsne.args$perplexity
        if (is.null(tsne.perplexity)) {
            tsne.perplexity <- formals(runTsne)$perplexity 
        }
        k.choices$runTsne <- tsnePerplexityToNeighbors(tsne.perplexity)
    }

    if (!is.null(clusterGraph.args)) {
        snn.k <- buildSnnGraph.args$num.neighbors
        if (is.null(snn.k)) {
            snn.k <- formals(buildSnnGraph)$num.neighbors
        }
        k.choices$clusterGraph <- snn.k
    }

    if (length(k.choices) == 0L) {
        output <- list()
        names(output) <- character(0)
        return(output)
    }

    if (is.matrix(x)) {
        index <- buildIndex(x, BNPARAM=BNPARAM, transposed=TRUE)
    } else {
        index <- x
    }
    nn.res <- list()

    if (collapse.search) {
        all.res <- findKNN(
            index,
            k=max(unlist(k.choices)),
            num.threads=num.threads,
            get.index="transposed",
            get.distance="transposed"
        )
        for (n in names(k.choices)) {
            curk <- k.choices[[n]]
            curi <- all.res$index
            curd <- all.res$distance
            if (curk != nrow(curi)) {
                curi <- curi[seq_len(curk),,drop=FALSE]
                curd <- curd[seq_len(curk),,drop=FALSE]
            }
            nn.res[[n]] <- list(index=curi, distance=curd)
        }

    } else {
        precomputed <- list()
        for (n in names(k.choices)) {
            curk <- k.choices[[n]]
            key <- as.character(curk) 
            if (key %in% names(precomputed)) {
                curres <- precomputed[[key]]
            } else {
                curres <- findKNN(
                    index,
                    k=curk,
                    num.threads=num.threads,
                    get.index="transposed",
                    get.distance="transposed"
                )
                precomputed[[key]] <- curres
            }
            nn.res[[n]] <- curres
        }
    }

    jobs <- list()
    if ("runUmap" %in% names(nn.res)) {
        jobs$runUmap <- list(runUmap=c(list(x=nn.res$runUmap), runUmap.args))
    }
    if ("runTsne" %in% names(nn.res)) {
        if (!("num.threads" %in% names(runTsne.args))) {
            runTsne.args$num.threads <- max(1L, num.threads - length(nn.res) + 1L)
        }
        jobs$runTsne <- list(runTsne=c(list(x=nn.res$runTsne), runTsne.args))
    }
    if ("clusterGraph" %in% names(nn.res)) {
        jobs$clusterGraph <- list(
            buildSnnGraph=c(list(x=nn.res$clusterGraph), buildSnnGraph.args),
            clusterGraph=clusterGraph.args,
            return.graph=return.graph
        )
    }

    # Don't attempt to be smart and use a FORK cluster on Unix-likes to squeeze
    # out more efficiency, as forking does not play nice with OpenMP; see
    # https://stackoverflow.com/questions/49049388/understanding-openmp-shortcomings-regarding-fork.
    # Just keep the PSOCK default.
    clust <- makeCluster(min(length(jobs), num.threads))
    on.exit(stopCluster(clust), add=TRUE, after=FALSE)
    output <- parLapply(clust, jobs, fun=.internal_executor)

    if (return.graph && "clusterGraph" %in% names(output)) {
        tmp <- output$clusterGraph
        output$clusterGraph <- NULL
        c(output, tmp)
    } else {
        output
    }
}

.internal_executor <- function(job) {
    if ("runUmap" %in% names(job)) {
        do.call(runUmap, job$runUmap)
    } else if ("runTsne" %in% names(job)) {
        do.call(runTsne, job$runTsne)
    } else {
        snn.graph <- do.call(buildSnnGraph, job$buildSnnGraph)
        clust <- do.call(clusterGraph, c(list(snn.graph), job$clusterGraph))

        if (job$return.graph) {
            list(buildSnnGraph=snn.graph, clusterGraph=clust)
        } else {
            return(clust)
        }
    }
}
