cmake <- "cmake"
options <- character(0)

if (.Platform$OS.type != "windows") {
    options <- c(options, "-DCMAKE_BUILD_TYPE=Release")
}

# Choosing the right compiler.
r.self <- file.path(R.home("bin"), "R")

c_compiler <- sub(" .*", "", system2(r.self, c("CMD", "CONFIG", "CC"), stdout=TRUE))
if (Sys.which(c_compiler) != "" || file.exists(c_compiler)) {
    options <- c(options, paste0("-DCMAKE_C_COMPILER=", c_compiler))
}

cxx_compiler <- sub(" .*", "", system2(r.self, c("CMD", "CONFIG", "CXX"), stdout=TRUE))
if (Sys.which(cxx_compiler) != "" || file.exists(cxx_compiler)) {
    options <- c(options, paste0("-DCMAKE_CXX_COMPILER=", cxx_compiler))
}

make <- sub(" .*", "", system2(r.self, c("CMD", "CONFIG", "MAKE"), stdout=TRUE))
if (Sys.which(make) != "" || file.exists(make)) {
    options <- c(options, paste0("-DCMAKE_MAKE_PROGRAM=", make))
}

# Setting up BLAS and LAPACK.
X <- sessionInfo()
blas_path <- X$BLAS
lapack_path <- X$LAPACK

if (file.exists(blas_path)) {
    options <- c(options,
        sprintf('-DBLAS_LIBRARIES="%s"', blas_path),
        "-DIGRAPH_USE_INTERNAL_BLAS=0"
    )
}

if (file.exists(lapack_path)) {
    options <- c(options,
        sprintf('-DLAPACK_LIBRARIES="%s"', lapack_path),
        "-DIGRAPH_USE_INTERNAL_LAPACK=0"
    )
}

# Forcing vendored builds of everything else.
options <- c(options,
    "-DIGRAPH_USE_INTERNAL_ARPACK=1",
    "-DIGRAPH_USE_INTERNAL_GLPK=1",
    "-DIGRAPH_USE_INTERNAL_GMP=1",
    "-DIGRAPH_USE_INTERNAL_PLFIT=1"
)

# Skipping the optional dependencies.
options <- c(options,
    "-DIGRAPH_BISON_SUPPORT=0",
    "-DIGRAPH_FLEX_SUPPORT=0"
)

# Downloading the file:
version <- "0.10.13"

build_path <- paste0("build_", version)
if (!file.exists(build_path)) {
    source_path <- paste0("source_", version)
    if (!file.exists(source_path)) {
        igraph_path <- "igraph-0.10.13.tar.gz"
        if (!file.exists(igraph_path)) {
            if (download.file("https://github.com/igraph/igraph/releases/download/0.10.13/igraph-0.10.13.tar.gz", igraph_path)) {
                stop("failed to download the igraph sources")
            }
        }
        untar(igraph_path, exdir=source_path)
    }

    options <- c(options, paste0("-DCMAKE_INSTALL_PREFIX=", file.path(normalizePath(dirname(getwd())), "src", "_deps")))
    system2(cmake, c("-S", file.path(source_path, paste0("igraph-", version)), "-B", build_path, options))
}

system2(cmake, c("--build", build_path))
system2(cmake, c("--install", build_path))
