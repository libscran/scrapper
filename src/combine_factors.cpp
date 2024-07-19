//#include "config.h"

#include "Rcpp.h"

#include <vector>
#include <algorithm>

#include "scran_aggregate/scran_aggregate.hpp"

static Rcpp::List convert_to_index_list(const std::vector<std::vector<int> >& levels) {
    size_t nfac = levels.size();
    Rcpp::List combos(nfac);
    for (size_t f = 0; f < nfac; ++f) {
        const auto& current = levels[f];
        combos[f] = Rcpp::IntegerVector(current.begin(), current.end());
    }
    return combos;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List combine_factors(Rcpp::List factors, bool keep_unused, Rcpp::IntegerVector nlevels) {
    size_t nfac = factors.size();
    if (nfac == 0) {
        throw std::runtime_error("'factors' must have length greater than zero");
    }

    std::vector<Rcpp::IntegerVector> ibuffers;
    ibuffers.reserve(nfac);
    for (size_t f = 0; f < nfac; ++f) {
        ibuffers.emplace_back(factors[f]);
    }

    size_t ngenes = ibuffers.front().size();
    for (size_t f = 1; f < nfac; ++f) {
        if (ibuffers[f].size() != ngenes) {
            throw std::runtime_error("all elements of 'factors' must have the same length");
        }
    }

    Rcpp::IntegerVector oindices;
    Rcpp::List olevels;

    if (keep_unused) {
        if (nlevels.size() != nfac) {
            throw std::runtime_error("'nlevels' and 'factors' must have the same length");
        }

        if (nfac > 1) {
            std::vector<std::pair<const int*, int> > buffers;
            buffers.reserve(nfac);
            for (size_t f = 0; f < nfac; ++f) {
                buffers.emplace_back(ibuffers[f].begin(), nlevels[f]);
            }
            oindices = Rcpp::IntegerVector(ngenes);
            auto res = scran_aggregate::combine_factors_unused(ngenes, buffers, static_cast<int*>(oindices.begin()));
            olevels = convert_to_index_list(res);
        } else {
            oindices = ibuffers[0];
            Rcpp::IntegerVector seq(nlevels[0]);
            std::iota(seq.begin(), seq.end(), 0);
            olevels = Rcpp::List::create(seq);
        }

    } else {
        if (nfac > 1) {
            std::vector<const int*> buffers;
            buffers.reserve(nfac);
            for (size_t f = 0; f < nfac; ++f) {
                buffers.emplace_back(ibuffers[f].begin());
            }
            oindices = Rcpp::IntegerVector(ngenes);
            auto res = scran_aggregate::combine_factors(ngenes, buffers, static_cast<int*>(oindices.begin()));
            olevels = convert_to_index_list(res);
        } else {
            oindices = Rcpp::IntegerVector(ngenes);
            auto res = scran_aggregate::clean_factor(ngenes, static_cast<const int*>(ibuffers[0].begin()), static_cast<int*>(oindices.begin()));
            olevels = Rcpp::List::create(Rcpp::IntegerVector(res.begin(), res.end()));
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("levels") = olevels,
        Rcpp::Named("index") = oindices
    );
}
