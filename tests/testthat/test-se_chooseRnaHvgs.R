# library(testthat); library(scrapple); source("test-se_chooseRnaHvgs.R")

set.seed(12341234)
library(SummarizedExperiment)
mean <- runif(100, 0, 5)
mat <- matrix(rpois(1000, mean * 10), ncol=10)
se <- SummarizedExperiment(list(counts=mat))
se <- normalizeRnaCounts.se(se)

test_that("chooseRnaHvgs.se works as expected", {
    out <- chooseRnaHvgs.se(se, top=20)
    expect_false(is.null(rowData(out)$means))
    expect_true(any(rowData(out)$hvg))
    expect_false(all(rowData(out)$hvg))

    # Works with a prefix.
    out <- chooseRnaHvgs.se(se, top=2000, output.prefix="VAR.")
    expect_type(rowData(out)$VAR.means, "double")
    expect_type(rowData(out)$VAR.variances, "double")
    expect_true(all(rowData(out)$VAR.hvgs))

    # Works with blocking.
    block <- rep(LETTERS[1:2], each=5)
    out <- chooseRnaHvgs.se(se, top=20, block=block)
    expect_null(rowData(out)$per.block)

    out <- chooseRnaHvgs.se(se, top=20, block=block, include.per.block=TRUE)
    expect_s4_class(rowData(out)$per.block, "DFrame")
    expect_type(rowData(out)$per.block$A$means, "double")
    expect_type(rowData(out)$per.block$B$residuals, "double")
})
