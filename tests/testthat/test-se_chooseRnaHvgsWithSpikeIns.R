# library(testthat); library(scrapper); source("test-se_chooseRnaHvgsWithSpikeIns.R")

set.seed(42424242)
library(SingleCellExperiment)
mat <- matrix(rpois(1000, 5), ncol=10)
sce <- SingleCellExperiment(list(counts=mat))
emat <- matrix(rpois(200, 5), ncol=10)
altExp(sce, "ERCC") <- SingleCellExperiment(list(counts=emat))
sce <- normalizeRnaCountsWithSpikeIns.se(sce, "ERCC")

test_that("chooseRnaHvgs.se works as expected", {
    out <- chooseRnaHvgsWithSpikeIns.se(sce, "ERCC", top=20)

    expect_type(rowData(out)$means, "double")
    expect_type(rowData(out)$residuals, "double")
    expect_equal(rowData(out)$fitted + rowData(out)$residuals, rowData(out)$variances)
    expect_true(any(rowData(out)$hvg))
    expect_false(all(rowData(out)$hvg))

    spike.rd <- rowData(altExp(out, "ERCC"))
    expect_type(spike.rd$means, "double")
    expect_type(spike.rd$residuals, "double")
    expect_null(spike.rd$hvg)

    # Same results with other methods of specifying the alternative experiment and assays. 
    alt <- chooseRnaHvgsWithSpikeIns.se(sce, c(ERCC=2L), top=20, assay.type=2)
    expect_identical(out, alt)

    # Works with a prefix.
    pre <- chooseRnaHvgsWithSpikeIns.se(sce, "ERCC", top=2000, output.prefix="VAR.")
    expect_type(rowData(pre)$VAR.means, "double")
    expect_type(rowData(altExp(pre))$VAR.variances, "double")
    expect_true(all(rowData(pre)$VAR.hvgs))
})
