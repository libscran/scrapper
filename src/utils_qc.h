#ifndef UTILS_QC_H
#define UTILS_QC_H

#include "Rcpp.h"

#include <vector>
#include <stdexcept>
#include <cstdint>

inline void cast_subset_vectors(size_t ngenes, const Rcpp::List& subsets, std::vector<Rcpp::LogicalVector>& store, std::vector<const int*>& ptrs) {
    size_t nsub = subsets.size();
    store.reserve(nsub);
    ptrs.reserve(nsub);
    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        store.emplace_back(*sIt);
        const auto& current = store.back();
        if (ngenes != static_cast<size_t>(current.size())) {
            throw std::runtime_error("each entry of 'subsets' should have the same length as 'nrow(x)'");
        }
        ptrs.emplace_back(current.begin());
    }
}

inline void prepare_subset_metrics(size_t ncells, size_t nsubsets, std::vector<Rcpp::NumericVector>& store, std::vector<double*>& ptrs) {
    store.reserve(nsubsets);
    ptrs.reserve(nsubsets);
    for (size_t s = 0; s < nsubsets; ++s) {
        store.emplace_back(ncells);
        ptrs.emplace_back(store.back().begin());
    }
}

inline void check_subset_metrics(size_t ncells, const Rcpp::List& input, std::vector<Rcpp::NumericVector>& store) {
    size_t nsubs = input.size();
    store.reserve(nsubs);
    for (size_t s = 0; s < nsubs; ++s) {
        store.emplace_back(input[s]);
        if (static_cast<size_t>(store.back().size()) != ncells) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }
    }
}

inline Rcpp::List create_subset_filters(const std::vector<std::vector<double> >& filters) {
    size_t nsubs = filters.size();
    Rcpp::List subs(nsubs);
    for (size_t s = 0; s < nsubs; ++s) {
        const auto& current = filters[s];
        subs[s] = Rcpp::NumericVector(current.begin(), current.end());
    }
    return subs;
}

inline void copy_filters_blocked(size_t nblocks, const Rcpp::NumericVector& input, std::vector<double>& store) {
    if (static_cast<size_t>(input.size()) != nblocks) {
        throw std::runtime_error("each array of thresholds in 'filters' should have length equal to the number of blocks");
    }
    store.insert(store.end(), input.begin(), input.end());
}

inline void copy_subset_filters_blocked(size_t nsubs, size_t nblocks, const Rcpp::List& subsets, std::vector<std::vector<double> >& store) {
    if (static_cast<size_t>(subsets.size()) != nsubs) {
        throw std::runtime_error("'filters.subsets' should have the same length as the number of subsets in 'metrics'");
    }
    store.resize(nsubs);
    for (size_t s = 0; s < nsubs; ++s) {
        Rcpp::NumericVector cursub(subsets[s]);
        copy_filters_blocked(nblocks, cursub, store[s]);
    }
}

inline double parse_filter_unblocked(const Rcpp::NumericVector& val, const char* msg) {
    if (val.size() != 1) {
        throw std::runtime_error("'" + std::string(msg) + "' should contain a single threshold");
    }
    return val[0]; 
}

inline void copy_subset_filters_unblocked(size_t nsubs, const Rcpp::NumericVector& subsets, std::vector<double>& store) {
    if (static_cast<size_t>(subsets.size()) != nsubs) {
        throw std::runtime_error("'filters.subsets' should have the same length as the number of subsets in 'metrics'");
    }
    store.insert(store.end(), subsets.begin(), subsets.end());
}

#endif
