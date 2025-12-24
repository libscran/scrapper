cache <- new.env()
cache$rna <- list()
cache$adt <- list()
cache$crispr <- list()

#' Get a test scRNA-seq dataset
#'
#' Get a test single-cell dataset with varying levels of processing.
#' This uses caching to avoid recomputation.
#'
#' @param at String specifying the level of processing.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing a dataset at the specified level of processing.
#' 
#' @details
#' For \code{getTestRnaData}, this is a scRNA-seq dataset of the mouse brain,
#' where the main experiment contains RNA counts and the alternative experiments contain ERCC and repeat element counts.
#' This is obtained with \code{\link[scRNAseq]{fetchDataset}("zeisel-brain-2015", "2023-12-14")}.
#'
#' For \code{getTestAdtData}, this is a CITE-seq dataset of human PBMCs,
#' where the main experiment contains RNA counts and the alternative experiment contains ADT counts.
#' This is obtained with \code{\link[scRNAseq]{fetchDataset}("kotliarov-pbmc-2020", "2024-04-18")}.
#'
#' For \code{getTestCrisprData}, this is a Perturb-seq dataset of a pancreatic beta cell line,
#' where the main experiment contains RNA counts and the alternative experiment contains CRISPR guide counts.
#' This is obtained with \code{\link[scRNAseq]{fetchDataset}("cao-pancreas-2025", "2025-10-10", "rqc")}.
#'
#' @author Aaron Lun
#' @examples
#' getTestRnaData.se()
#' getTestAdtData.se()
#' getTestCrisprData.se()
#'
#' @seealso
#' \code{\link[scRNAseq]{fetchDataset}}, used to obtain each dataset.
#' 
#' @export
#' @name getTestData.se
getTestRnaData.se <- function(at = c("start", "qc", "norm", "hvg", "pca", "cluster")) {
    at <- match.arg(at)

    if (!("start" %in% names(cache$rna))) {
        cache$rna$start <- scRNAseq::fetchDataset("zeisel-brain-2015", "2023-12-14", realize.assays=TRUE)
    }
    sce <- cache$rna$start
    if (at == "start") {
        return(sce)
    }

    if (!("qc" %in% names(cache$rna))) {
        sce <- quickRnaQc.se(sce, subsets=list(mito=startsWith(rownames(sce), "mt-")), altexp.proportions="ERCC")
        sce <- sce[,sce$keep]
        cache$rna$qc <- sce
    }
    sce <- cache$rna$qc
    if (at == "qc") {
        return(sce)
    }

    if (!("norm" %in% names(cache$rna))) {
        sce <- normalizeRnaCounts.se(sce)
        cache$rna$norm <- sce
    }
    sce <- cache$rna$norm
    if (at == "norm") {
        return(sce)
    }

    if (!("hvg" %in% names(cache$rna))) {
        sce <- chooseRnaHvgs.se(sce, more.var.args=list(use.min.width=TRUE))
        cache$rna$hvg <- sce
    }
    sce <- cache$rna$hvg
    if (at == "hvg") {
        return(sce)
    }

    if (!("pca" %in% names(cache$rna))) {
        sce <- runPca.se(sce, features=rowData(sce)$hvg)
        cache$rna$pca <- sce
    }
    sce <- cache$rna$pca
    if (at == "pca") {
        return(sce)
    }

    if (!("cluster" %in% names(cache$rna))) {
        sce <- clusterGraph.se(sce)
        cache$rna$cluster <- sce
    }
    sce <- cache$rna$cluster
    if (at == "cluster") {
        return(sce)
    }
}

#' @export
#' @rdname getTestData.se
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom SingleCellExperiment altExp altExp<-
getTestAdtData.se <- function(at = c("start", "qc", "norm", "hvg", "pca")) {
    at <- match.arg(at)

    if (!("start" %in% names(cache$adt))) {
        raw.sce <- scRNAseq::fetchDataset("kotliarov-pbmc-2020", "2024-04-18")
        raw.sce <- raw.sce[,1:5000] # Cutting it down a bit for speed.
        assay(raw.sce) <- as(assay(raw.sce), "dgCMatrix")
        assay(altExp(raw.sce)) <- as(assay(altExp(raw.sce)), "dgCMatrix")
        cache$adt$start <- raw.sce
    }
    sce <- cache$adt$start
    if (at == "start") {
        return(sce)
    }

    if (!("qc" %in% names(cache$adt))) {
        sce <- quickRnaQc.se(sce, subsets=list(mito=startsWith(rownames(sce), "MT-")))
        altExp(sce, "ADT") <- quickAdtQc.se(altExp(sce, "ADT"), subsets=list(igg=rowData(altExp(sce, "ADT"))$isotype))
        sce <- sce[,sce$keep & altExp(sce, "ADT")$keep]
        cache$adt$qc <- sce
    }
    sce <- cache$adt$qc
    if (at == "qc") {
        return(sce)
    }

    if (!("norm" %in% names(cache$adt))) {
        sce <- normalizeRnaCounts.se(sce)
        altExp(sce, "ADT") <- normalizeAdtCounts.se(altExp(sce, "ADT"))
        cache$adt$norm <- sce
    }
    sce <- cache$adt$norm
    if (at == "norm") {
        return(sce)
    }

    if (!("hvg" %in% names(cache$adt))) {
        sce <- chooseRnaHvgs.se(sce)
        cache$adt$hvg <- sce
    }
    sce <- cache$adt$hvg
    if (at == "hvg") {
        return(sce)
    }

    if (!("pca" %in% names(cache$adt))) {
        sce <- runPca.se(sce, features=rowData(sce)$hvg)
        altExp(sce, "ADT") <- runPca.se(altExp(sce, "ADT"), features=NULL)
        cache$adt$pca <- sce
    }
    sce <- cache$adt$pca
    if (at == "pca") {
        return(sce)
    }
}

#' @export
#' @rdname getTestData.se
#' @importFrom methods as
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom SingleCellExperiment altExp altExp<-
getTestCrisprData.se <- function(at = c("start", "qc")) {
    at <- match.arg(at)

    if (!("start" %in% names(cache$crispr))) {
        raw.sce <- scRNAseq::fetchDataset("cao-pancreas-2025", "2025-10-10", "rqc")
        raw.sce <- raw.sce[,1:5000] # Cutting it down a bit for speed.
        assay(raw.sce) <- as(assay(raw.sce), "dgCMatrix")
        assay(altExp(raw.sce)) <- as(assay(altExp(raw.sce)), "dgCMatrix")
        cache$crispr$start <- raw.sce
    }
    sce <- cache$crispr$start
    if (at == "start") {
        return(cache$crispr[[at]])
    }

    if (!("qc" %in% names(cache$crispr))) {
        sce <- quickRnaQc.se(sce, subsets=list(mito=startsWith(rownames(sce), "MT-")))
        altExp(sce, "CRISPR Guide Capture") <- quickCrisprQc.se(altExp(sce, "CRISPR Guide Capture"))
        sce <- sce[,sce$keep & altExp(sce, "CRISPR Guide Capture")$keep]
        cache$crispr$qc <- sce
    }
    sce <- cache$crispr$qc
    if (at == "qc") {
        return(sce)
    }
}
