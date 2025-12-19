#ifndef UTILS_QC_H
#define UTILS_QC_H

#include "config.h"

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <unordered_map>

#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

template<typename Ngenes_>
void cast_subset_vectors(Ngenes_ ngenes, const Rcpp::List& subsets, std::vector<Rcpp::LogicalVector>& store, std::vector<const int*>& ptrs) {
    const auto nsub = subsets.size();
    store.reserve(nsub);
    ptrs.reserve(nsub);

    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        store.emplace_back(*sIt);
        const auto& current = store.back();
        if (!sanisizer::is_equal(ngenes, current.size())) {
            throw std::runtime_error("each entry of 'subsets' should have the same length as 'nrow(x)'");
        }
        ptrs.emplace_back(current.begin());
    }
}

template<typename Ncells_, typename Nsubsets_>
void prepare_subset_metrics(Ncells_ ncells, Nsubsets_ nsubsets, std::vector<Rcpp::NumericVector>& store, std::vector<double*>& ptrs) {
    store.reserve(nsubsets);
    ptrs.reserve(nsubsets);
    sanisizer::as_size_type<Rcpp::NumericVector>(ncells);
    for (I<decltype(nsubsets)> s = 0; s < nsubsets; ++s) {
        store.emplace_back(ncells);
        ptrs.emplace_back(store.back().begin());
    }
}

template<typename Ncells_>
void check_subset_metrics(Ncells_ ncells, const Rcpp::List& input, std::vector<Rcpp::NumericVector>& store) {
    const auto nsubs = input.size();
    store.reserve(nsubs);
    for (I<decltype(nsubs)> s = 0; s < nsubs; ++s) {
        store.emplace_back(input[s]);
        if (!sanisizer::is_equal(store.back().size(), ncells)) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }
    }
}

inline Rcpp::List create_subset_filters(const std::vector<std::vector<double> >& filters) {
    const auto nsubs = filters.size();
    Rcpp::List subs(nsubs);
    for (I<decltype(nsubs)> s = 0; s < nsubs; ++s) {
        const auto& current = filters[s];
        subs[s] = Rcpp::NumericVector(current.begin(), current.end());
    }
    return subs;
}

template<typename Nblocks_>
void copy_filters_blocked(Nblocks_ nblocks, const Rcpp::NumericVector& input, std::vector<double>& store) {
    if (!sanisizer::is_equal(input.size(), nblocks)) {
        throw std::runtime_error("each array of thresholds in 'filters' should have length equal to the number of blocks");
    }
    store.insert(store.end(), input.begin(), input.end());
}

template<typename Nsubs_, typename Nblocks_>
void copy_subset_filters_blocked(Nsubs_ nsubs, Nblocks_ nblocks, const Rcpp::List& subsets, std::vector<std::vector<double> >& store) {
    if (!sanisizer::is_equal(subsets.size(), nsubs)) {
        throw std::runtime_error("'filters.subsets' should have the same length as the number of subsets in 'metrics'");
    }
    sanisizer::resize(store, nsubs);
    for (I<decltype(nsubs)> s = 0; s < nsubs; ++s) {
        Rcpp::NumericVector cursub(subsets[s]);
        copy_filters_blocked(nblocks, cursub, store[s]);
    }
}

inline void check_names(const Rcpp::List& contents, const std::vector<std::string>& names) {
    if (!contents.hasAttribute("names")) {
        throw std::runtime_error("list is unnamed");
    }
    Rcpp::CharacterVector all_names(contents.attr("names"));

    std::unordered_set<std::string> known_names;
    for (auto x : all_names) {
        known_names.insert(Rcpp::as<std::string>(x));
    }

    for (const auto& n : names) {
        if (known_names.find(n) == known_names.end()) {
            throw std::runtime_error("expected list to contain element with name '" + n + "'");
        }
    }
}

inline double parse_filter_unblocked(const Rcpp::NumericVector& val, const char* msg) {
    if (val.size() != 1) {
        throw std::runtime_error("'" + std::string(msg) + "' should contain a single threshold");
    }
    return val[0]; 
}

template<typename Nsubs_>
void copy_subset_filters_unblocked(Nsubs_ nsubs, const Rcpp::NumericVector& subsets, std::vector<double>& store) {
    if (!sanisizer::is_equal(subsets.size(), nsubs)) {
        throw std::runtime_error("'filters.subsets' should have the same length as the number of subsets in 'metrics'");
    }
    store.insert(store.end(), subsets.begin(), subsets.end());
}

#endif
