\name{scrapper News}
\title{scrapper News}
\encoding{UTF-8}

\section{Version 1.2}{\itemize{
\item Added the \code{aggregateAcrossGenes()} function, to compute an aggregate expression value for gene sets.

\item Added \code{compute.cohens.d=}, \code{compute.delta.mean=} and \code{compute.delta.detected=} options to \code{scoreMarkers()}.

\item Support \code{top=Inf} in \code{chooseHighlyVariableGenes()}.
Also added the \code{bound=} argument to set a hard upper/lower bound.

\item Bugfix for correct filtering with \code{block=} in the various \code{filter*QcMetrics()} functions when not all blocking levels are present.

\item Bugfix to \code{clusterKmeans()} to respect the user-supplied \code{seed=} in relevant initialization methods.

\item Added a \code{return.graph=} option to return the SNN graph from \code{runAllNeighborSteps()}.

\item Added a \code{testEnrichment()} function for quick and dirty gene set enrichment testing.
}}
