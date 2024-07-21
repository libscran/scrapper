# library(testthat); library(scrapper); source("test-fitVarianceTrend.R")

test_that("fitVarianceTrend works as expected", {
    x <- runif(1000)
    y <- 2^rnorm(1000)
    out <- fitVarianceTrend(x, y)
    expect_equal(y - out$fitted, out$residuals)

    out2 <- fitVarianceTrend(x, y, num.threads=2)
    expect_identical(out, out2)

    # Responds to the various options.
    out2 <- fitVarianceTrend(x, y, use.min.width=TRUE, min.width=0.5)
    expect_false(identical(out, out2))

    out2 <- fitVarianceTrend(x, y, transform=FALSE)
    expect_false(identical(out, out2))

    out2 <- fitVarianceTrend(x, y, span=0.5)
    expect_false(identical(out, out2))
})
