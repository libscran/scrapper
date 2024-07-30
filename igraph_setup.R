cmake <- "cmake"
options <- character(0)

options <- c(options, 
     "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
     "-DIGRAPH_WARNINGS_AS_ERRORS=OFF"
)

if (.Platform$OS.type != "windows") {
    options <- c(options, "-DCMAKE_BUILD_TYPE=Release")
}

if (Sys.info()[["sysname"]] == "Darwin") {
     options <- c(options, "-DCMAKE_OSX_DEPLOYMENT_TARGET=\"\"") # avoiding hard-coding of exact OSX version.
}

# Choosing the right compiler.
r.self <- file.path(R.home("bin"), "R")

c_compiler <- sub(" .*", "", system2(r.self, c("CMD", "config", "CC"), stdout=TRUE))
if (Sys.which(c_compiler) != "" || file.exists(c_compiler)) {
    options <- c(options, paste0("-DCMAKE_C_COMPILER=", c_compiler))
}

cxx_compiler <- sub(" .*", "", system2(r.self, c("CMD", "config", "CXX"), stdout=TRUE))
if (Sys.which(cxx_compiler) != "" || file.exists(cxx_compiler)) {
    options <- c(options, paste0("-DCMAKE_CXX_COMPILER=", cxx_compiler))
}

make <- sub(" .*", "", system2(r.self, c("CMD", "config", "MAKE"), stdout=TRUE))
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

install_path <- file.path("src", "_deps")

if (!file.exists(install_path)) {
    tmp_dir <- "igraph_temp"
    dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)
    build_path <- file.path(tmp_dir, paste0("build_", version))

    if (!file.exists(build_path)) {
        source_path <- file.path(tmp_dir, paste0("source_", version))
        if (!file.exists(source_path)) {
            igraph_path <- file.path(tmp_dir, paste0("igraph-", version, ".tar.gz"))

            if (!file.exists(igraph_path)) {
                if (download.file(paste0("https://github.com/igraph/igraph/releases/download/", version, "/igraph-", version, ".tar.gz"), igraph_path)) {
                    stop("failed to download the igraph sources")
                }
            }
            untar(igraph_path, exdir=source_path)
        }

        options <- c(options, paste0("-DCMAKE_INSTALL_PREFIX=", install_path))
        system2(cmake, c("-S", file.path(source_path, paste0("igraph-", version)), "-B", build_path, options))
    }

    if (.Platform$OS.type != "windows") {
        system2(cmake, c("--build", build_path))
    } else {
        system2(cmake, c("--build", build_path, "--config", "Release"))
    }

    system2(cmake, c("--install", build_path), stdout=FALSE)
}
