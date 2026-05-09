#include "config.h"

#include <cstddef>

#include "qdtsne/qdtsne.hpp"

#include "utils_neighbors.h"
#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP run_tsne(
    Rcpp::IntegerMatrix nnidx,
    Rcpp::NumericMatrix nndist,
    Rcpp::RObject perplexity,
    Rcpp::RObject theta,
    Rcpp::RObject early_exaggeration_iterations,
    Rcpp::RObject exaggeration_factor,
    Rcpp::RObject momentum_switch_iterations,
    Rcpp::RObject start_momentum,
    Rcpp::RObject final_momentum,
    Rcpp::RObject eta,
    Rcpp::RObject max_depth,
    Rcpp::RObject leaf_approximation,
    Rcpp::RObject max_iterations,
    double seed,
    Rcpp::RObject num_threads
) {
    qdtsne::Options opt;
    set_number(perplexity, opt.perplexity, "perplexity");
    opt.precomputed_perplexity_policy = qdtsne::PrecomputedPerplexityPolicy::CHECK;
    set_number(theta, opt.theta, "theta");
    set_integer(early_exaggeration_iterations, opt.early_exaggeration_iterations, "early.exaggeration.iterations");
    set_number(exaggeration_factor, opt.exaggeration_factor, "exaggeration.factor");
    set_integer(momentum_switch_iterations, opt.momentum_switch_iterations, "momentum.switch.iterations");
    set_number(start_momentum, opt.start_momentum, "start.momentum");
    set_number(final_momentum, opt.final_momentum, "final.momentum");
    set_number(eta, opt.eta, "eta");
    set_bool(leaf_approximation, opt.leaf_approximation, "leaf.approximation");
    set_integer(max_depth, opt.max_depth, "max.depth");
    set_integer(max_iterations, opt.max_iterations, "max.iterations");
    set_integer(num_threads, opt.num_threads, "num.threads");

    auto neighbors = unpack_neighbors<int, double>(nnidx, nndist);
    const auto nobs = neighbors.size();
    auto status = qdtsne::initialize<2>(std::move(neighbors), opt);

    auto output = create_matrix<Rcpp::NumericMatrix>(2, nobs);
    auto optr = static_cast<double*>(output.begin());
    qdtsne::initialize_random<2>(optr, sanisizer::cast<std::size_t>(nobs), sanisizer::from_float<unsigned long long>(seed));
    status.run(optr);

    return output;
}

//[[Rcpp::export(rng=false)]]
int perplexity_to_neighbors(double p) {
    return qdtsne::perplexity_to_k(p);
}

//[[Rcpp::export(rng=false)]]
Rcpp::List run_tsne_defaults() {
    Rcpp::List output;
    qdtsne::Options opt;
    output["perplexity"] = opt.perplexity;
    output["theta"] = opt.theta;
    output["early.exaggeration.iterations"] = opt.early_exaggeration_iterations;
    output["exaggeration.factor"] = opt.exaggeration_factor;
    output["momentum.switch.iterations"] = opt.momentum_switch_iterations;
    output["start.momentum"] = opt.start_momentum;
    output["final.momentum"] = opt.final_momentum;
    output["eta"] = opt.eta;
    output["leaf.approximation"] = opt.leaf_approximation;
    output["max.depth"] = opt.max_depth;
    output["max.iterations"] = opt.max_iterations;
    output["num.threads"] = opt.num_threads;
    return output;
}
