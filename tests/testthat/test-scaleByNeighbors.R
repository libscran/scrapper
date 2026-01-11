# library(testthat); library(scrapper); source("test-scaleByNeighbors.R")

set.seed(2000)
pcs <- list(
    gene = matrix(rnorm(10000), ncol=200),
    protein = matrix(rnorm(1000, sd=3), ncol=200),
    guide = matrix(rnorm(2000, sd=5), ncol=200)
)

manual_knn_distance <- function(x, k = 20) {
    value <- numeric(ncol(x))
    for (i in seq_len(ncol(x))) {
        chosen <- x[,i]
        distances <- sqrt(colSums((x[,-i] - chosen)^2))
        value[i] <- sort(distances)[k]
    } 
    median(value)
}

test_that("scaleByNeighbors works as expected", {
    out <- scaleByNeighbors(pcs, BNPARAM=BiocNeighbors::VptreeParam())
    expect_identical(names(out$scaling), names(pcs))

    references <- lapply(pcs, manual_knn_distance)
    expect_equal(out$scaling, references[[1]]/unlist(references))

    expect_identical(ncol(out$combined), ncol(pcs[[1]]))
    expect_identical(nrow(out$combined), sum(vapply(pcs, nrow, 0L)))

    expect_error(scaleByNeighbors(SummarizedExperiment::SummarizedExperiment()), "not supported")
})

test_that("scaleByNeighbors works as expected with blocking", {
    # Same results with duplicated blocks.
    {
        shuffle <- sample(200)
        block <- rep(LETTERS[1:3], each=200)[shuffle]
        pcs3 <- lapply(pcs, function(x) {
            cbind(x, x, x)[,shuffle,drop=FALSE]
        })

        ref <- scaleByNeighbors(pcs, BNPARAM=BiocNeighbors::VptreeParam())
        out <- scaleByNeighbors(pcs3, block=block, BNPARAM=BiocNeighbors::VptreeParam())
        expect_identical(out$scaling, ref$scaling)
    }

    # Manual comparison using a randomly sampled block.
    {
        set.seed(20001)
        block <- sample(LETTERS[1:3], 200, replace=TRUE)
        out <- scaleByNeighbors(pcs, block=block, block.weight.policy="equal", BNPARAM=BiocNeighbors::VptreeParam())

        references <- lapply(pcs, function(x) {
            cur.d <- numeric()
            for (b in LETTERS[1:3]) {
                cur.d <- append(cur.d, manual_knn_distance(x[,b == block,drop=FALSE]))
            }
            mean(cur.d)
        })

        references <- 1/unlist(references)
        normalized <- references / (references[1] / out$scaling[1])
        expect_equal(normalized, out$scaling)
    }
})
