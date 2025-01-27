//#include "config.h"

#include "Rcpp.h"
#include "knncolle/knncolle.hpp"
#include "qdtsne/qdtsne.hpp"

#include "utils_neighbors.h"

//[[Rcpp::export(rng=false)]]
SEXP run_tsne(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist, double perplexity, int leaf_approx, int max_depth, int max_iter, int seed, int num_threads) {
    qdtsne::Options opt;
    opt.perplexity = perplexity;
    opt.infer_perplexity = false; // rely on the perplexity supplied by the user.
    opt.leaf_approximation = leaf_approx;
    opt.max_depth = max_depth;
    opt.max_iterations = max_iter;
    opt.num_threads = num_threads;

    auto neighbors = unpack_neighbors<int, double>(nnidx, nndist);
    size_t nobs = neighbors.size();
    auto status = qdtsne::initialize<2>(std::move(neighbors), opt);

    Rcpp::NumericMatrix output(2, nobs);
    auto optr = static_cast<double*>(output.begin());
    qdtsne::initialize_random<2>(optr, nobs, seed);
    status.run(optr);

    return output;
}

//[[Rcpp::export(rng=false)]]
int perplexity_to_neighbors(double p) {
    return qdtsne::perplexity_to_k(p);
}
