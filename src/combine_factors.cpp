#include "config.h"

#include "Rcpp.h"

#include <vector>

#include "scran_aggregate/scran_aggregate.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::List combine_factors(Rcpp::List factors) {
    size_t nfac = factors.size();

    std::vector<Rcpp::IntegerVector> ibuffers;
    std::vector<const int*> buffers;
    ibuffers.reserve(nfac);
    buffers.reserve(nfac);
    for (size_t f = 0; f < nfac; ++f) {
        ibuffers.emplace_back(factors[f]);
        buffers.push_back(ibuffers.begin());
    }

    if (ibuffers.size()) {
        throw std::runtime_error("at least one factor must be supplied");
    }
    size_t ngenes = ibuffers.front().size();
    for (size_t f = 1; f < nfac; ++f) {
        if (ibuffers[f].size() != ngenes) {
            throw std::runtime_error("all factors must have the same length");
        }
    }

    Rcpp::IntegerVector output(ngenes);
    auto res = scran_aggregate::combine_factors(ngenes, buffers, static_cast<int*>(output.begin()));

    Rcpp::List combos(nfac);
    for (size_t f = 0; f < nfac; ++f) {
        const auto& curfac = ref.factors[f];
        combos[f] = Rcpp::IntegerVector(curfac.begin(), curfac.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("levels") = combos,
        Rcpp::Named("index") = output
    );
}
