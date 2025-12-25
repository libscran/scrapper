# library(testthat); library(scrapper); source("test-se_runPca.R")

library(SingleCellExperiment)
set.seed(8888)
sce <- SingleCellExperiment(list(logcounts = matrix(rnorm(100000), 500, 200)))

test_that("runPca.se works as expected", {
    out <- runPca.se(sce, features=1:50) 
    expect_true("PCA" %in% reducedDimNames(out))
    expect_identical(nrow(metadata(out)$PCA$rotation), 50L)
    expect_identical(length(metadata(out)$PCA$variance.explained), 25L)

    # Works with features set to NULL.
    null <- runPca.se(sce[1:50,], features=NULL) 
    expect_identical(reducedDim(null), reducedDim(out))
    expect_identical(metadata(null)$PCA$rotation, metadata(out)$PCA$rotation)

    no.meta <- runPca.se(sce, features=1:10, meta.name=NULL)
    expect_null(metadata(no.meta)$PCA)
})

test_that("runPca.se works with a basic SE", {
    se <- SummarizedExperiment(list(logcounts=assay(sce)))
    out <- runPca.se(se, features=5:100) 
    expect_s4_class(out, "SingleCellExperiment")
    expect_identical(reducedDimNames(out), "PCA")
})
