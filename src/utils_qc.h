#ifndef UTILS_QC_H
#define UTILS_QC_H

#include "config.h"

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <unordered_set>
#include <optional>

#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

template<typename Vector_, int sexptype_, class Message_>
Vector_ safe_cast(const Rcpp::RObject& raw_val, Message_ msg) {
    if (raw_val.sexp_type() != sexptype_) {
        std::string type;
        switch (sexptype_) {
            case INTSXP:
                type = "integer vector";
            case REALSXP:
                type = "double-precision vector";
            case LGLSXP:
                type = "logical vector";
            default:
                type = "list";
        }
        throw std::runtime_error(msg() + " should be a " + type);
    }
    return Vector_(raw_val);
}

template<typename Vector_, int sexptype_>
Vector_ parse_metrics(const Rcpp::RObject& raw_val, const char* field) {
    return safe_cast<Vector_, sexptype_>(raw_val, [&]() -> std::string { return "'metrics$" + std::string(field) + "'"; });
}

template<typename Ngenes_>
void cast_subset_vectors(Ngenes_ ngenes, const Rcpp::List& subsets, std::vector<Rcpp::LogicalVector>& store, std::vector<const int*>& ptrs) {
    const auto nsub = subsets.size();
    store.reserve(nsub);
    ptrs.reserve(nsub);

    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        store.emplace_back(safe_cast<Rcpp::LogicalVector, LGLSXP>(*sIt, []() -> std::string { return "each entry of 'subsets'"; }));
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
void check_subset_metrics(Ncells_ ncells, const Rcpp::RObject& raw_input, std::vector<Rcpp::NumericVector>& store) {
    auto input = safe_cast<Rcpp::List, VECSXP>(raw_input, []() -> std::string { return "'metrics$subset'"; });
    const auto nsubs = input.size();
    store.reserve(nsubs);

    for (I<decltype(nsubs)> s = 0; s < nsubs; ++s) {
        store.emplace_back(safe_cast<Rcpp::NumericVector, REALSXP>(input[s], []() -> std::string { return "each entry of 'metrics$subsets'"; }));
        if (!sanisizer::is_equal(store.back().size(), ncells)) {
            throw std::runtime_error("length of each 'metrics$subsets' vector should be equal to the number of cells");
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

template<typename Nblocks_, class Message_>
void copy_filters_blocked_raw(std::optional<Nblocks_>& nblocks, const Rcpp::RObject& raw_input, Message_ msg, std::vector<double>& store) {
    auto input = safe_cast<Rcpp::NumericVector, REALSXP>(raw_input, msg);
    if (nblocks.has_value()){
        if (!sanisizer::is_equal(input.size(), *nblocks)) {
            throw std::runtime_error(msg() + " should have length equal to the number of blocks");
        }
    } else {
        nblocks = sanisizer::cast<Nblocks_>(input.size());
    }
    store.insert(store.end(), input.begin(), input.end());
}

template<typename Nblocks_>
void copy_filters_blocked(std::optional<Nblocks_>& nblocks, const Rcpp::RObject& raw_input, const char* field, std::vector<double>& store) {
    copy_filters_blocked_raw(nblocks, raw_input, [&]() -> std::string { return "'thresholds$" + std::string(field) + "'"; }, store);
}

template<typename Nsubs_, typename Nblocks_>
void copy_subset_filters_blocked(Nsubs_ nsubs, Nblocks_ nblocks, const Rcpp::RObject& raw_subsets, std::vector<std::vector<double> >& store) {
    auto subsets = safe_cast<Rcpp::List, VECSXP>(raw_subsets, []() -> std::string { return "'thresholds$subset'"; });
    if (!sanisizer::is_equal(subsets.size(), nsubs)) {
        throw std::runtime_error("'thresholds$subsets' should have the same length as 'metrics$subsets'");
    }
    sanisizer::resize(store, nsubs);
    for (I<decltype(nsubs)> s = 0; s < nsubs; ++s) {
        copy_filters_blocked_raw(nblocks, subsets[s], []() -> std::string { return "each entry of 'thresholds$subsets'"; }, store[s]);
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

inline double parse_filter_unblocked(const Rcpp::RObject& raw_val, const char* msg) {
    auto val = safe_cast<Rcpp::NumericVector, REALSXP>(raw_val, [&]() -> std::string { return "'" + std::string(msg) + "'"; });
    if (val.size() != 1) {
        throw std::runtime_error("'" + std::string(msg) + "' should contain a single threshold");
    }
    return val[0]; 
}

template<typename Nsubs_>
void copy_subset_filters_unblocked(Nsubs_ nsubs, const Rcpp::RObject& raw_subsets, std::vector<double>& store) {
    auto subsets = safe_cast<Rcpp::NumericVector, REALSXP>(raw_subsets, [&]() -> std::string { return "'thresholds$subsets'"; });
    if (!sanisizer::is_equal(subsets.size(), nsubs)) {
        throw std::runtime_error("'thresholds$subsets' should have the same length as the number of subsets in 'metrics'");
    }
    store.insert(store.end(), subsets.begin(), subsets.end());
}

#endif
