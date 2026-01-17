#include "config.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "factorize/factorize.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

static Rcpp::List convert_to_index_list(const std::vector<std::vector<int> >& levels) {
    const auto nfac = levels.size();
    auto combos = sanisizer::create<Rcpp::List>(nfac);
    for (I<decltype(nfac)> f = 0; f < nfac; ++f) {
        const auto& current = levels[f];
        combos[f] = Rcpp::IntegerVector(current.begin(), current.end());
    }
    return combos;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List combine_factors(Rcpp::List factors, bool keep_unused, Rcpp::IntegerVector nlevels) {
    const auto nfac = factors.size();
    if (nfac == 0) {
        throw std::runtime_error("'factors' must have length greater than zero");
    }

    std::vector<Rcpp::IntegerVector> ibuffers;
    ibuffers.reserve(nfac);
    for (I<decltype(nfac)> f = 0; f < nfac; ++f) {
        ibuffers.emplace_back(factors[f]);
    }

    const auto ngenes = ibuffers.front().size();
    for (I<decltype(nfac)> f = 1; f < nfac; ++f) {
        if (!sanisizer::is_equal(ibuffers[f].size(), ngenes)) {
            throw std::runtime_error("all elements of 'factors' must have the same length");
        }
    }

    Rcpp::IntegerVector oindices;
    Rcpp::List olevels;

    if (keep_unused) {
        if (!sanisizer::is_equal(nlevels.size(), nfac)) {
            throw std::runtime_error("'nlevels' and 'factors' must have the same length");
        }

        std::vector<std::pair<const int*, int> > buffers;
        buffers.reserve(nfac);
        for (I<decltype(nfac)> f = 0; f < nfac; ++f) {
            buffers.emplace_back(ibuffers[f].begin(), nlevels[f]);
        }

        oindices = sanisizer::create<Rcpp::IntegerVector>(ngenes);
        auto res = factorize::combine_to_factor_unused(ngenes, buffers, static_cast<int*>(oindices.begin()));
        olevels = convert_to_index_list(res);

    } else {
        std::vector<const int*> buffers;
        buffers.reserve(nfac);
        for (I<decltype(nfac)> f = 0; f < nfac; ++f) {
            buffers.emplace_back(ibuffers[f].begin());
        }

        oindices = Rcpp::IntegerVector(ngenes);
        auto res = factorize::combine_to_factor(ngenes, buffers, static_cast<int*>(oindices.begin()));
        olevels = convert_to_index_list(res);
    }

    return Rcpp::List::create(
        Rcpp::Named("levels") = olevels,
        Rcpp::Named("index") = oindices
    );
}
