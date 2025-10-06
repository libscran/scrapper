//#include "config.h"

#include "Rcpp.h"
#include "knncolle/knncolle.hpp"
#include "qdtsne/qdtsne.hpp"

#include "utils_neighbors.h"

//[[Rcpp::export(rng=false)]]
SEXP run_tsne(
    Rcpp::IntegerMatrix nnidx,
    Rcpp::NumericMatrix nndist,
    double perplexity,
    double theta,
    int early_exaggeration_iterations,
    double exaggeration_factor,
    int momentum_switch_iterations,
    double start_momentum,
    double final_momentum,
    double eta,
    int max_depth,
    int leaf_approx,
    int max_iter,
    int seed,
    int num_threads)
{
    qdtsne::Options opt;
    opt.perplexity = perplexity;
    opt.infer_perplexity = false; // rely on the perplexity supplied by the user.
    opt.theta = theta;
    opt.early_exaggeration_iterations = early_exaggeration_iterations;
    opt.exaggeration_factor = exaggeration_factor;
    opt.momentum_switch_iterations = momentum_switch_iterations;
    opt.start_momentum = start_momentum;
    opt.final_momentum = final_momentum;
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
