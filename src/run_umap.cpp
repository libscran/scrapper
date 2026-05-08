#include "config.h"

#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "umappp/umappp.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_neighbors.h"
#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP run_umap(
    Rcpp::IntegerMatrix nnidx,
    Rcpp::NumericMatrix nndist,
    int num_dim,
    Rcpp::RObject local_connectivity,
    Rcpp::RObject bandwidth,
    Rcpp::RObject mix_ratio,
    Rcpp::RObject spread,
    Rcpp::RObject min_dist,
    Rcpp::RObject a,
    Rcpp::RObject b,
    Rcpp::RObject repulsion_strength,
    Rcpp::RObject initialize_method,
    Rcpp::Nullable<Rcpp::NumericMatrix> initial_coordinates,
    Rcpp::RObject initialize_random_on_spectral_fail,
    Rcpp::RObject initialize_spectral_scale,
    Rcpp::RObject initialize_spectral_jitter,
    Rcpp::RObject initialize_spectral_jitter_sd,
    Rcpp::RObject initialize_random_scale,
    Rcpp::RObject initialize_seed,
    Rcpp::RObject num_epochs,
    Rcpp::RObject learning_rate,
    Rcpp::RObject negative_sample_rate,
    Rcpp::RObject optimize_seed,
    Rcpp::RObject num_threads,
    Rcpp::RObject num_threads_optimize 
) {
    auto neighbors = unpack_neighbors<int, float>(nnidx, nndist);
    const auto nobs = neighbors.size();

    umappp::Options opt;
    set_number(local_connectivity, opt.local_connectivity, "local.connectivity");
    set_number(bandwidth, opt.bandwidth, "bandwidth");
    set_number(mix_ratio, opt.mix_ratio, "mix.ratio");
    set_number(spread, opt.spread, "spread");
    set_number(min_dist, opt.min_dist, "min.dist");
    set_optional_number(a, opt.a, "a");
    set_optional_number(b, opt.b, "b");
    set_number(repulsion_strength, opt.repulsion_strength, "repulsion.strength");

    if (!initialize_method.isNULL()) {
        const std::string im = parse_single_string(initialize_method, "initialize.method");
        if (im == "spectral") {
            opt.initialize_method = umappp::InitializeMethod::SPECTRAL;
        } else if (im == "random") {
            opt.initialize_method = umappp::InitializeMethod::RANDOM;
        } else if (im == "none") {
            opt.initialize_method = umappp::InitializeMethod::NONE;
        } else {
            throw std::runtime_error("unknown value for 'initialize_method'");
        }
    }

    set_bool(initialize_random_on_spectral_fail, opt.initialize_random_on_spectral_fail, "initialize.random.on.spectral.fail");
    set_number(initialize_spectral_scale, opt.initialize_spectral_scale, "initialize.spectral.scale");
    set_bool(initialize_spectral_jitter, opt.initialize_spectral_jitter, "initialize.spectral.jitter");
    set_number(initialize_spectral_jitter_sd, opt.initialize_spectral_jitter_sd, "initialize.spectral.jitter.sd");
    set_number(initialize_random_scale, opt.initialize_random_scale, "initialize.random.scale");
    set_integer(initialize_seed, opt.initialize_seed, "initialize.seed");

    std::vector<float> embedding(sanisizer::product<typename std::vector<float>::size_type>(num_dim, nobs));
    if (!initial_coordinates.isNull()) {
        Rcpp::NumericMatrix realized(initial_coordinates);
        std::copy(realized.begin(), realized.end(), embedding.data());
    } else if (opt.initialize_method == umappp::InitializeMethod::NONE || !opt.initialize_random_on_spectral_fail) {
        throw std::runtime_error("expected initial coordinates to be supplied");
    }

    set_optional_integer(num_epochs, opt.num_epochs, "num.epochs");
    set_number(learning_rate, opt.learning_rate, "learning.rate");
    set_number(negative_sample_rate, opt.negative_sample_rate, "negative.sample.rate");
    set_integer(optimize_seed, opt.optimize_seed, "optimize.seed");
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_integer(num_threads_optimize, opt.num_threads_optimize, "num.threads.optimize");

    auto status = umappp::initialize(std::move(neighbors), sanisizer::cast<std::size_t>(num_dim), embedding.data(), opt);
    status.run(embedding.data());

    auto output = create_matrix<Rcpp::NumericMatrix>(num_dim, nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List run_umap_defaults() {
    Rcpp::List output;

    umappp::Options opt;
    output["num.neighbors"] = opt.num_neighbors;
    output["local.connectivity"] = opt.local_connectivity;
    output["bandwidth"] = opt.bandwidth;
    output["mix.ratio"] = opt.mix_ratio;
    output["spread"] = opt.spread;
    output["min.dist"] = opt.min_dist;

    if (!opt.a.has_value()) {
        output["a"] = R_NilValue;
    } else {
        throw std::runtime_error("unexpected a default for runUmap");
    }

    if (!opt.b.has_value()) {
        output["b"] = R_NilValue;
    } else {
        throw std::runtime_error("unexpected b default for runUmap");
    }

    output["repulsion.strength"] = opt.repulsion_strength;

    if (opt.initialize_method == umappp::InitializeMethod::SPECTRAL) {
        output["initialize.method"] = "spectral";
    } else {
        throw std::runtime_error("unexpected initialize_method for runUmap");
    }

    output["initialize.random.on.spectral.fail"] = opt.initialize_random_on_spectral_fail;
    output["initialize.spectral.scale"] = opt.initialize_spectral_scale;
    output["initialize.spectral.jitter"] = opt.initialize_spectral_jitter;
    output["initialize.spectral.jitter.sd"] = opt.initialize_spectral_jitter_sd;
    output["initialize.random.scale"] = opt.initialize_random_scale;
    output["initialize.seed"] = opt.initialize_seed;

    if (!opt.num_epochs.has_value()) {
        output["num.epochs"] = R_NilValue;
    } else {
        throw std::runtime_error("unexpected num_epochs for runUmap");
    }

    output["learning.rate"] = opt.learning_rate;
    output["negative.sample.rate"] = opt.negative_sample_rate;
    output["optimize.seed"] = opt.optimize_seed;
    output["num.threads"] = opt.num_threads;
    output["num.threads.optimize"] = opt.num_threads_optimize;

    return output;
}
