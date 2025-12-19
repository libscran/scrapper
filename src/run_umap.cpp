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
    double local_connectivity,
    double bandwidth,
    double mix_ratio,
    double spread,
    double min_dist,
    Rcpp::Nullable<Rcpp::NumericVector> a,
    Rcpp::Nullable<Rcpp::NumericVector> b,
    double repulsion_strength,
    std::string initialize_method,
    Rcpp::Nullable<Rcpp::NumericMatrix> initial_coordinates,
    bool initialize_random_on_spectral_fail,
    double initialize_spectral_scale,
    bool initialize_spectral_jitter,
    double initialize_spectral_jitter_sd,
    double initialize_random_scale,
    double initialize_seed,
    Rcpp::Nullable<Rcpp::IntegerVector> num_epochs,
    double learning_rate,
    double negative_sample_rate,
    double optimize_seed,
    int num_threads,
    bool parallel_optimization)
{
    auto neighbors = unpack_neighbors<int, float>(nnidx, nndist);
    const auto nobs = neighbors.size();

    umappp::Options opt;
    opt.local_connectivity = local_connectivity;
    opt.bandwidth = bandwidth;
    opt.mix_ratio = mix_ratio;
    opt.spread = spread;
    opt.min_dist = min_dist;

    if (!a.isNull()) {
        Rcpp::NumericVector realized(a);
        if (realized.size() != 1) {
            throw std::runtime_error("'a' should be a numeric scalar or NULL");
        }
        opt.a = realized[0];
    }

    if (!b.isNull()) {
        Rcpp::NumericVector realized(b);
        if (realized.size() != 1) {
            throw std::runtime_error("'b' should be a numeric scalar or NULL");
        }
        opt.b = realized[0];
    }

    opt.repulsion_strength = repulsion_strength;

    if (initialize_method == "spectral") {
        opt.initialize_method = umappp::InitializeMethod::SPECTRAL;
    } else if (initialize_method == "random") {
        opt.initialize_method = umappp::InitializeMethod::RANDOM;
    } else if (initialize_method == "none") {
        opt.initialize_method = umappp::InitializeMethod::NONE;
    } else {
        throw std::runtime_error("unknown value for 'initialize_method'");
    }

    std::vector<float> embedding(sanisizer::product<typename std::vector<float>::size_type>(num_dim, nobs));
    if (!initial_coordinates.isNull()) {
        Rcpp::NumericMatrix realized(initial_coordinates);
        std::copy(realized.begin(), realized.end(), embedding.data());
    } else if (initialize_method == "none" || !initialize_random_on_spectral_fail) {
        throw std::runtime_error("expected initial coordinates to be supplied");
    }

    opt.initialize_random_on_spectral_fail = initialize_random_on_spectral_fail;
    opt.initialize_spectral_scale = initialize_spectral_scale;
    opt.initialize_spectral_jitter = initialize_spectral_jitter;
    opt.initialize_spectral_jitter_sd = initialize_spectral_jitter_sd;
    opt.initialize_random_scale = initialize_random_scale;
    opt.initialize_seed = sanisizer::from_float<decltype(opt.initialize_seed)>(initialize_seed);

    if (!num_epochs.isNull()) {
        Rcpp::IntegerVector realized(num_epochs);
        if (realized.size() != 1) {
            throw std::runtime_error("'num_epochs' should be an integer scalar or NULL");
        }
        opt.num_epochs = realized[0];
    }

    opt.learning_rate = learning_rate;
    opt.negative_sample_rate = negative_sample_rate;
    opt.optimize_seed = sanisizer::from_float<decltype(opt.optimize_seed)>(optimize_seed);
    opt.num_threads = num_threads;
    opt.parallel_optimization = parallel_optimization;

    auto status = umappp::initialize(std::move(neighbors), sanisizer::cast<std::size_t>(num_dim), embedding.data(), opt);
    status.run(embedding.data());

    auto output = create_matrix<Rcpp::NumericMatrix>(num_dim, nobs);
    std::copy(embedding.begin(), embedding.end(), output.begin());
    return output;
}
