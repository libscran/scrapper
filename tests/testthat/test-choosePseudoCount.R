# library(testthat); library(scrapper); source("test-choosePseudoCount.R")

set.seed(9999)
test_that("choosePseudoCount works as expected", {
    sf <- 2^rnorm(100, sd=2)
    out <- choosePseudoCount(sf)
    expect_gt(out, 0)
    expect_lt(out, choosePseudoCount(sf, quantile=0.01))
    expect_lt(out, choosePseudoCount(sf, max.bias=0.5))
})

test_that("defaults work correctly", {
    def <- choosePseudoCountDefaults()
    expect_true(all(names(def) %in% names(formals(choosePseudoCount))))
})
