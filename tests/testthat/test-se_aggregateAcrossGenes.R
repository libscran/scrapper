# library(testthat); library(scrapper); source("test-se_aggregateAcrossGenes.R")

library(SummarizedExperiment)
mat <- matrix(runif(1000), ncol=10)
se <- SummarizedExperiment(list(logcounts=mat))
se$stuff <- sample(LETTERS[1:5], ncol(se), replace=TRUE)
se$whee <- sample(c(TRUE, FALSE), ncol(se), replace=TRUE)

sets <- list(foo=sample(nrow(mat), 10), bar=sample(nrow(mat), 20), stuff=sample(nrow(mat), 30))

test_that("aggregateAcrossGenes.se works as expected", {
    out <- aggregateAcrossGenes.se(se, sets)
    expect_identical(assayNames(out), "logcounts")
    expect_identical(rownames(out), names(sets))
    expect_identical(colData(out), colData(se))

    # Checking the names.
    out2 <- aggregateAcrossGenes.se(se, sets, assay.type=1)
    expect_identical(assayNames(out2), "aggregated")
    out2 <- aggregateAcrossGenes.se(se, sets, assay.type=1, output.name="FOO")
    expect_identical(assayNames(out2), "FOO")

    # Works with empty inputs.
    zeroed <- aggregateAcrossGenes.se(se, sets[0])
    expect_identical(nrow(zeroed), 0L)
})

test_that("aggregateAcrossGenes.se works with Lists", {
    out <- aggregateAcrossGenes.se(se, sets)

    lsets <- List(sets)
    lout <- aggregateAcrossGenes.se(se, lsets)
    expect_identical(lout, out)

    mcols(lsets) <- DataFrame(description=c("fooish", "barish", "stuffish"))
    lout2 <- aggregateAcrossGenes.se(se, lsets)
    expect_identical(rowData(lout2), mcols(lsets))

    # Cleans up DataFrames.
    wsets <- lapply(sets, function(x) DataFrame(index=x, weight=runif(length(x))))
    wout <- aggregateAcrossGenes.se(se, wsets)
    wref <- aggregateAcrossGenes.se(se, lapply(wsets, as.list))
    expect_identical(wout, wref)
})
