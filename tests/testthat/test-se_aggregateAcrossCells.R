# library(testthat); library(scrapper); source("test-se_aggregateAcrossCells.R")

library(SingleCellExperiment)
mat <- matrix(rpois(1000, 10), ncol=100)
se <- SummarizedExperiment(list(counts=mat))
se$stuff <- sample(LETTERS[1:5], ncol(se), replace=TRUE)
se$whee <- sample(c(TRUE, FALSE), ncol(se), replace=TRUE)
rowData(se)$foo <- runif(nrow(se))

test_that("aggregateAcrossCells.se works as expected", {
    out <- aggregateAcrossCells.se(se, colData(se)[,"stuff",drop=FALSE])
    expect_identical(ncol(out), length(unique(se$stuff)))
    expect_identical(rowData(out), rowData(se))
    expect_identical(out$factor.stuff, out$stuff)
    expect_identical(sort(out$factor.stuff), sort(unique(se$stuff)))
    expect_identical(out$counts, as.integer(table(se$stuff)[out$stuff]))
    expect_identical(metadata(out)$aggregated$index, match(se$stuff, out$stuff))

    # Works with direct specification of the factor.
    out2 <- aggregateAcrossCells.se(se, se$stuff)
    expect_identical(out2$factor.1, out2$stuff)
    expect_identical(assay(out2), assay(out))

    # Works if we provide a list.
    out3 <- aggregateAcrossCells.se(se, list(foo=se$stuff), output.prefix=NULL)
    expect_identical(out3$foo, out$stuff)
    expect_identical(assay(out3), assay(out))

    # Works with multiple factors.
    mult <- aggregateAcrossCells.se(se, colData(se)[,c("stuff", "whee"),drop=FALSE])
    expect_identical(mult$factor.stuff, mult$stuff)
    expect_identical(mult$factor.whee, mult$whee)
    counters <- S4Vectors::match(colData(se)[,c("stuff", "whee"),drop=FALSE], colData(mult)[,c("stuff", "whee"),drop=FALSE])
    expect_identical(as.integer(table(counters)), mult$counts)
    expect_identical(metadata(mult)$aggregated$index, counters)

    # Works if we disable the outputs.
    out4 <- aggregateAcrossCells.se(se, se$stuff, meta.name=NULL, counts.name=NULL, include.coldata=FALSE)
    expect_null(out4$counts)
    expect_null(metadata(out4)$aggregated)
    expect_null(out4$stuff)
})

test_that("aggregateAcrossCells.se works with alternative experiments", {
    sce <- as(se, "SingleCellExperiment")
    is.even <- 1:10 %% 2 == 0
    altExp(sce, "GuP") <- se[is.even,]
    sce$random <- paste0("FOO_", sce$stuff)

    out <- aggregateAcrossCells.se(sce, colData(se)[,c("stuff", "whee"),drop=FALSE])
    expect_false(is(out, "SingleCellExperiment"))

    out <- aggregateAcrossCells.se(sce, colData(se)[,c("stuff", "whee"),drop=FALSE], altexps="GuP")
    expect_s4_class(out, "SingleCellExperiment")
    expect_identical(assay(out)[is.even,], assay(altExp(out)))
    expect_identical(rowData(out)[is.even,,drop=FALSE], rowData(altExp(out)))
    expect_identical(out$stuff, altExp(out)$stuff)
    expect_identical(out$random, paste0("FOO_", out$stuff))
    expect_null(altExp(out)$factor.stuff)
    expect_null(altExp(out)$factor.whee)
    expect_null(altExp(out)$counts)
    expect_null(metadata(altExp(out))$aggregated)
    expect_null(altExp(out)$random)

    # Copying restores the extra outputs.
    out <- aggregateAcrossCells.se(sce, colData(se)[,c("stuff", "whee"),drop=FALSE], altexps="GuP", copy.altexps=TRUE)
    expect_s4_class(out, "SingleCellExperiment")
    expect_identical(out$stuff, altExp(out)$factor.stuff)
    expect_identical(out$whee, altExp(out)$factor.whee)
    expect_identical(out$counts, altExp(out)$counts)
    expect_identical(metadata(out)$aggregated, metadata(altExp(out))$aggregated)
})

test_that("aggregateColData works as expected", {
    df <- DataFrame(stuff=rep(LETTERS[1:5], each=4), whee=rep(c(TRUE, FALSE), 10))

    adf <- aggregateColData(df, as.integer(factor(df$stuff)), 5)
    expect_identical(adf$stuff, LETTERS[1:5])
    expect_identical(adf$whee, rep(NA, 5))

    adf <- aggregateColData(df, as.integer(factor(df$whee)), 2)
    expect_identical(adf$stuff, rep(NA_character_, 2))
    expect_identical(adf$whee, c(FALSE, TRUE))

    # Still works as a factor.
    df$stuff <- factor(df$stuff)
    adf <- aggregateColData(df, as.integer(factor(df$stuff)), 5)
    expect_identical(adf$stuff, factor(LETTERS[1:5]))

    # Handles non-atomic, non-factor columns.
    df$foo <- DataFrame(bar=runif(nrow(df)))
    adf <- aggregateColData(df, as.integer(factor(df$whee)), 2)
    expect_null(adf$foo)
    adf <- aggregateColData(df, as.integer(factor(df$whee)), 2, only.atomic=FALSE)
    expect_identical(adf$foo, rep(NA, 2))

    zeroed <- aggregateColData(df[0,], integer(0), 0)
    expect_true(is.factor(zeroed$stuff))
    expect_type(zeroed$whee, "logical")
})
