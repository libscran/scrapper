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
\item Output now contains the number of rows, row names and group IDs for easier interpretation by downstream functions.
}
\item Added arguments to \code{summarizeEffects()} to disable calculation of certain summary statistics, e.g., by setting \code{compute.summary.median=FALSE}.
Also support calculation of arbitrary quantiles as described for \code{scoreMarkers()}.
\item \code{suggest*QcThresholds()} functions now return a \code{block.ids} field that preserves the type of \code{block=}.
This enables easier matching in the corresponding \code{filter*QcMetrics()} function.
\item Multiple changes to \code{modelGeneVariances()}:
\itemize{
\item Allow the use of quantiles to average statistics across blocks via the new \code{block.average.policy=} and \code{block.quantile=} arguments.
\item Return a \code{block.ids} field that preserves the type of \code{block=} for easier matching.
}
\item Multiple changes to \code{runPca()}:
\itemize{
\item Added a \code{subset=} argument that performs a PCA on a subset of features but still reports the full rotation matrix for all features. 
\item Return a \code{block.ids} field that preserves the type of \code{block=} for easier matching.
}
\item \code{clusterKmeans()} now emits a warning upon convergence failure.
This can be disabled by setting \code{warn=FALSE}.
\item Simplified the automatic naming performed by \code{combineFactors()} when \code{factors} is unnamed.
\item Added SummarizedExperiment-compatible wrappers for all functions.
These can be identified by their \code{*.se} suffix, e.g., \code{aggregateAcrossCells.se()}.
Some functions execute multiple steps for convenience, e.g., \code{quickRnaQc.se()} runs all of the RNA-related QC functions.
\item Added \code{analyze.se()} to replace the \code{analyze()} function.
The former is more convenient it stores most outputs in a SingleCellExperiment and formats the marker tables correctly.
The latter is now deprecated.
\item \code{compute*QcMetrics()} functions now return a DataFrame for more convenient inspection and manipulation.
\item All functions that previously returned a base data.frame will now return a DataFrame from \pkg{S4Vectors}.
This is a little prettier and we already load \pkg{S4Vectors} anyway.
\item Removed checks on \code{num.neighbors=} when pre-computed neighbor search results are supplied in \code{runUmap()}, \code{buildSnnGraph()} and \code{subsampleByNeighbors()}.
Any differences between \code{num.neighbors=} and the number of neighbors in the pre-computed results are now ignored.
\item Implemented the LogNormalizedMatrix class for delayed log-transformation.
Modified \code{normalizeCounts()} to return instances of this class when \code{log=TRUE} and \code{delayed=TRUE}.
This allows downstream C++ code to use a more efficient normalization wrapper via \code{initializeCpp()}.
Note that there are some very minor differences in the log-normalized values when computed in R via \code{extract_array()} and those in C++ via \pkg{tatami} extraction.
\item Added the \code{countGroupsByBlock()} function to count the number of cells within each combination of group and block.
This is primarily intended for diagnosing batch effects.
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
