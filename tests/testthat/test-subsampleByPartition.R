# library(testthat); library(scrapper); source("test-subsampleByPartition.R")

set.seed(2000)

test_that("subsampleByPartition works as expected", {
    part <- sample(LETTERS, 1000, replace=TRUE)
    keep <- subsampleByPartition(part, 100)
    expect_identical(length(keep), 100L)
    expect_true(all(keep >= 1))
    expect_true(all(keep <= 1000))
    expect_identical(sort(unique(part)), sort(unique(part[keep])))

    part <- c(LETTERS[1:10], sample(LETTERS[11:26], 990, replace=TRUE))
    keep <- subsampleByPartition(part, 100)
    expect_identical(sort(unique(part)), sort(unique(part[keep])))
    keep <- subsampleByPartition(part, 100, force.non.empty=FALSE)
    expect_false(all(LETTERS[1:10] %in% part[keep]))
})
