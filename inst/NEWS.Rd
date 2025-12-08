\name{scrapper News}
\title{scrapper News}
\encoding{UTF-8}

\section{Version 1.6.0}{\itemize{
\item Multiple changes to \code{clusterGraph()}:
\itemize{
\item Updates to use the latest version of the \pkg{igraph} library.
\item \code{method="leiden"} now supports the ER objective function in \code{leiden.objective=}.
\item Non-zero status are now directly reported as errors, \code{status} is removed from the output.
}
\item Multiple additions to \code{scoreMarkers()}:
\itemize{
\item Support use of quantiles to average statistics across blocks, via the new \code{block.average.policy=} and \code{block.quantile=} arguments.
\item Allow calculation of arbitrary quantiles as a summary statistic for the effect sizes of each group, via the new \code{compute.summary.quantiles=} argument.
}
\item Added arguments to \code{summarizeEffects()} to disable calculation of certain summary statistics, e.g., by setting \code{compute.summary.median=FALSE}.
Also support calculation of arbitrary quantiles as described for \code{scoreMarkers()}.
\item Modified \code{modelGeneVariances()} to allow the use of quantiles to average statistics across blocks,
via the new \code{block.average.policy=} and \code{block.quantile=} arguments.
\item Added a \code{subset=} argument to \code{runPca()} to allow users to perform a PCA on a subset of features,
but still report the full rotation matrix for all features. 
}}

\section{Version 1.4.0}{\itemize{
\item Multiple changes to \code{scoreMarkers()}:
\itemize{
\item Setting \code{all.pairwise=} to an integer will now return the top markers by effect size for each pairwise comparison.
This is more efficient than creating the 3D array of effects and then ranking the genes in each comparison.
\item The min-rank summary is now limited to the top \code{min.rank.limit=} genes, for efficiency.
All other genes will have their reported min-rank set to a maximum placeholder.
\item Added \code{compute.group.mean=} and \code{compute.group.detected=} options to disable reporting of the group-wise means.
}

\item Multiple changes to \code{testEnrichment()}:
\itemize{
\item Setting \code{universe=NULL} will automatically define the universe as the union of all genes.
\item Return value is now a data frame with the number of overlapping genes and the size of each gene set.
}

\item Updated \code{correctMnn()} to the new algorithm used by the \pkg{mnncorrect} C++ library.
This aims to improve correction accuracy and reduce sensitivity to the choice of reference.

\item Bugfix to \code{runAllNeighborSteps()} to return an empty list when no step is requested.

\item Set \code{bound=0} by default in \code{chooseHighlyVariableGenes()}, for more convenient use with residuals.

\item Modified \code{scoreGeneSet()} to return \code{weights} as a data frame containing the row indices of the genes in the set.

\item Support character and logical vectors in \code{sets=} for \code{aggregateAcrossGenes()}.

\item Automatically remove empty clusters from the output of \code{clusterKmeans()}.

\item Supported blocking and related arguments for weighted averages of distances in \code{scaleByNeighbors()}.
}}

\section{Version 1.2.0}{\itemize{
\item Added the \code{aggregateAcrossGenes()} function, to compute an aggregate expression value for gene sets.

\item Added \code{compute.cohens.d=}, \code{compute.delta.mean=} and \code{compute.delta.detected=} options to \code{scoreMarkers()}.

\item Support \code{top=Inf} in \code{chooseHighlyVariableGenes()}.
Also added the \code{bound=} argument to set a hard upper/lower bound.

\item Bugfix for correct filtering with \code{block=} in the various \code{filter*QcMetrics()} functions when not all blocking levels are present.

\item Bugfix to \code{clusterKmeans()} to respect the user-supplied \code{seed=} in relevant initialization methods.

\item Added a \code{return.graph=} option to return the SNN graph from \code{runAllNeighborSteps()}.

\item Added a \code{testEnrichment()} function for quick and dirty gene set enrichment testing.

\item Modified \code{runPca()} so that it caps \code{number=} to the maximum number of available PCs.

\item Added an \code{analyze()} function that provides a one-click approach for analyzing single-cell data.

\item Added a \code{reportGroupMarkerStatistics()} function to combine all marker statistics for a single group into one data frame.

\item Switch \code{clusterGraph()} to use C++ wrappers around the igraph community detection algorithms via \pkg{Rigraphlib}.
This replaces the dependency on the igraph R package.

\item Added a \code{delayed=} option to avoid wrapping the matrix in a \code{DelayedArray} in \code{normalizeCounts()}.
}}
