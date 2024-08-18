//#include "config.h"

#include "Rcpp.h"
#include "umappp/umappp.hpp"

#include "utils_neighbors.h"

//[[Rcpp::export(rng=false)]]
SEXP run_umap(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist, int ndim, double min_dist, int seed, int num_epochs, int num_threads, bool parallel_optimization) {
    auto neighbors = unpack_neighbors<int, float>(nnidx, nndist);
    size_t nobs = neighbors.size();

    umappp::Options opt;
    opt.min_dist = min_dist;
    opt.seed = seed;
    opt.num_epochs = num_epochs;
    opt.num_threads = num_threads;
    opt.parallel_optimization = parallel_optimization;

    std::vector<float> embedding(ndim * nobs);
    auto status = umappp::initialize(std::move(neighbors), ndim, embedding.data(), opt);
    status.run();

    Rcpp::NumericMatrix output(ndim, nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}
