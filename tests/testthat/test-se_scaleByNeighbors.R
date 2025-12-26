# library(testthat); library(scrapper); source("test-se_scaleByNeighbors.R")

library(SingleCellExperiment)
set.seed(1111)
sce <- SingleCellExperiment(list(counts = matrix(0, 0, 200)))
reducedDim(sce, "PCA1") <- matrix(runif(1000), nrow=200)
reducedDim(sce, "PCA2") <- matrix(runif(1000), nrow=200)
altExp(sce, "ADT") <- sce
altExp(sce, "other") <- sce

test_that("scaleByNeighbors.se works as expected", {
    out <- scaleByNeighbors.se(
        sce,
        altexp.reddims=list(
            ADT=c("PCA2", "PCA1"),
            other=c("PCA1", "PCA2")
        ),
        main.reddims=c("PCA1", "PCA2")
    )

    ref <- scaleByNeighbors(
        list(
            t(reducedDim(sce, "PCA1")),
            t(reducedDim(sce, "PCA2")),
            t(reducedDim(altExp(sce, "ADT"), "PCA2")),
            t(reducedDim(altExp(sce, "ADT"), "PCA1")),
            t(reducedDim(altExp(sce, "other"), "PCA1")),
            t(reducedDim(altExp(sce, "other"), "PCA2"))
        )
    )
    expect_identical(reducedDim(out, "combined"), t(ref$combined))

    meta <- metadata(out)$combined
    expect_identical(names(meta$main.scaling), c("PCA1", "PCA2"))
    expect_identical(names(meta$altexp.scaling), c("ADT", "other"))
    expect_identical(names(meta$altexp.scaling$ADT), c("PCA2", "PCA1"))
    expect_identical(names(meta$altexp.scaling$other), c("PCA1", "PCA2"))

    expect_identical(
        c(
            meta$main.scaling[["PCA1"]],
            meta$main.scaling[["PCA2"]],
            meta$altexp.scaling[["ADT"]][["PCA2"]],
            meta$altexp.scaling[["ADT"]][["PCA1"]],
            meta$altexp.scaling[["other"]][["PCA1"]],
            meta$altexp.scaling[["other"]][["PCA2"]]
        ),
        ref$scaling
    )
})

test_that("scaleByNeighbors.se works with edge cases", {
    # Main experiment only.
    {
        out <- scaleByNeighbors.se(sce, altexp.reddims=list(), main.reddims=c("PCA2", "PCA1"))
        ref <- scaleByNeighbors(list(t(reducedDim(sce, "PCA2")), t(reducedDim(sce, "PCA1"))))

        expect_identical(reducedDim(out, "combined"), t(ref$combined))
        expect_identical(names(metadata(out)$combined$main.scaling), c("PCA2", "PCA1"))
        expect_identical(unname(metadata(out)$combined$main.scaling), ref$scaling)
    }

    # Alternative experiments only.
    {
        out <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT="PCA2", other="PCA2"), main.reddims=NULL)
        ref <- scaleByNeighbors(list(t(reducedDim(altExp(sce, "ADT"), "PCA2")), t(reducedDim(altExp(sce, "other"), "PCA2"))))
        expect_identical(reducedDim(out, "combined"), t(ref$combined))

        meta.alt <- metadata(out)$combined$altexp.scaling
        expect_identical(names(metadata(out)$combined$altexp.scaling), c("ADT", "other"))
        expect_identical(names(metadata(out)$combined$altexp.scaling$ADT), "PCA2")
        expect_identical(names(metadata(out)$combined$altexp.scaling$other), "PCA2")
        expect_identical(c(meta.alt[["ADT"]][["PCA2"]], meta.alt[["other"]][["PCA2"]]), ref$scaling)
    }

    # Skips the metadata if requested.
    out <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT="PCA2", other="PCA2"), main.reddims="PCA1", meta.name=NULL)
    expect_null(metadata(out)$combined)
})

test_that("scaleByNeighbors.se works deduplicates correctly", {
    out <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT=c("PCA1", "PCA1"), ADT=c("PCA2", "PCA2")), main.reddims=c("PCA2", "PCA2"))
    ref <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT=c("PCA1", "PCA1")), main.reddims="PCA2")
    expect_identical(reducedDim(out, "combined"), reducedDim(ref, "combined"))
    expect_identical(metadata(out)$combined, metadata(ref)$combined)
})
