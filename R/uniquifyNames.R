uniquifyNames <- function(old.names, new.names) {
    missing.name <- is.na(new.names)
    new.names[missing.name] <- old.names[missing.name]
    is.dup.name <- new.names %in% new.names[duplicated(new.names)]
    new.names[is.dup.name] <- paste0(new.names[is.dup.name], "_", old.names[is.dup.name])
    new.names
}
