//#include "config.h"

#include "Rcpp.h"

#include "nenesub/nenesub.hpp"
#include "tatami/tatami.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP subsample_by_neighbors(Rcpp::IntegerMatrix indices, Rcpp::NumericMatrix distances, int min_remaining) {
    const auto num_obs = distances.cols();
    const auto k = distances.rows();
    if (!sanisizer::is_equal(indices.rows(), k) || !sanisizer::is_equal(indices.cols(), num_obs)) {
        throw std::runtime_error("'indices' and 'distances' must have the same dimensions");
    }

    if (sanisizer::is_less_than(k, min_remaining)) {
        throw std::runtime_error("'min_remaining' should not be greater than the number of neighbors");
    }

    const int* iptr = indices.begin();
    const double* dptr = distances.begin();

    nenesub::Options opt;
    opt.min_remaining = min_remaining;
    std::vector<int> selected;

    nenesub::compute(
        num_obs,
        /* get_neighbors = */ [&](I<decltype(num_obs)> i) -> tatami::ArrayView<int> {
            return tatami::ArrayView<int>(iptr + sanisizer::product_unsafe<std::size_t>(k, i), k);
        },
        /* get_index = */ [](const tatami::ArrayView<int>& neighbors, I<decltype(num_obs)> i) -> int {
            return neighbors[i] - 1;
        },
        /* get_max_distance = */ [&](I<decltype(num_obs)> i) -> double {
            return dptr[sanisizer::nd_offset<std::size_t>(k - 1, k, i)];
        },
        opt, 
        selected
    );

    for (auto& s : selected) {
        ++s;
    }
    return Rcpp::IntegerVector(selected.begin(), selected.end());
}
