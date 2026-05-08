# library(testthat); library(scrapper); source("test-centerSpikeInFactors.R")

set.seed(8888)

endogenous <- runif(100)
spike.ercc <- runif(100)
spike.sirv <- runif(100)

test_that("centerSpikeInFactors works as expected", {
    centered <- centerSpikeInFactors(endogenous, list(ERCC = spike.ercc, SIRV = spike.sirv))
    expect_equal(centered$endogenous, endogenous / mean(endogenous))
    expect_equal(centered$spike.ins$ERCC, spike.ercc / mean(spike.ercc))
    expect_equal(centered$spike.ins$SIRV, spike.sirv / mean(spike.sirv))
})

test_that("centerSpikeInFactors works with blocking", {
    b <- sample(3, 100, replace=TRUE)

    centered <- centerSpikeInFactors(endogenous, list(ERCC = spike.ercc, SIRV = spike.sirv), block = b)
    e.scaling <- vapply(split(endogenous, b), mean, 0)
    expect_equal(centered$endogenous, endogenous / min(e.scaling))
    eout.scaling <- vapply(split(centered$endogenous, b), mean, 0)
    s1.scaling <- vapply(split(spike.ercc, b), mean, 0)
    expect_equal(centered$spike.ins$ERCC, spike.ercc / unname(s1.scaling[b] / eout.scaling[b]))
    s2.scaling <- vapply(split(spike.sirv, b), mean, 0)
    expect_equal(centered$spike.ins$SIRV, spike.sirv / unname(s2.scaling[b] / eout.scaling[b]))

    centered2 <- centerSpikeInFactors(endogenous, list(ERCC = spike.ercc, SIRV = spike.sirv), block = b, mode = "lowest")
    expect_identical(centered, centered2)

    centered <- centerSpikeInFactors(endogenous, list(ERCC = spike.ercc, SIRV = spike.sirv), block = b, mode = "per-block")
    e.scaling <- vapply(split(endogenous, b), mean, 0)
    expect_equal(centered$endogenous, endogenous / unname(e.scaling[b]))
    s1.scaling <- vapply(split(spike.ercc, b), mean, 0)
    expect_equal(centered$spike.ins$ERCC, spike.ercc / unname(s1.scaling[b]))
    s2.scaling <- vapply(split(spike.sirv, b), mean, 0)
    expect_equal(centered$spike.ins$SIRV, spike.sirv / unname(s2.scaling[b]))

    expect_error(centerSpikeInFactors(endogenous, list(ERCC = spike.ercc, SIRV = spike.sirv), block = b, mode = "WHEE"), "unknown")
})

test_that("defaults work correctly", {
    def <- centerSpikeInFactorsDefaults()
    expect_true(all(names(def) %in% names(formals(centerSpikeInFactors))))
})
