uniquifyVectorsByGroup <- function(x, grouping, drop.fail = TRUE) {
    by.group <- split(seq_along(grouping), grouping)
    output <- vector("list", length(x))

    for (i in seq_along(x)) {
        current <- x[[i]]

        uniq <- try({
            if (!is.null(dim(y))) {
                uvals <- lapply(x, function(y) {
                    got <- unique(current[y,,drop=FALSE])
                    if (dim(got)[1] != 1) {
                        NA
                    } else {
                        got
                    }
                })
                do.call(rbind, uvals)
            } else {
                uvals <- lapply(x, function(y) {
                    got <- unique(current[y])
                    if (length(got) != 1) {
                        NA
                    } else {
                        got
                    }
                })
                do.call(c, uvals)
            }
        }, silent=TRUE)

        if (is(uniq, 'try-error') || NROW(uniq) != length(x)) {
            if (!drop.fail) {
                output[[i]] <- rep(NA, length(x))
            }
        } else {
            output[[i]] <- uniq
        }
    }



}

