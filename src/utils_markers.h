#ifndef UTILS_MARKERS_H
#define UTILS_MARKERS_H

#include "config.h"

#include <vector>
#include <cstddef>
#include <string>
#include <optional>

#include "scran_markers/scran_markers.hpp"

#include "utils_other.h"

template<typename Ngroups_, typename Ngenes_, typename Nquantiles_>
inline void initialize_summary_buffers(
    const Ngroups_ num_groups,
    const Ngenes_ num_genes,
    std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs,
    const bool compute_min,
    std::vector<Rcpp::NumericVector>& min,
    const bool compute_mean,
    std::vector<Rcpp::NumericVector>& mean,
    const bool compute_median,
    std::vector<Rcpp::NumericVector>& median,
    const bool compute_max,
    std::vector<Rcpp::NumericVector>& max,
    const Nquantiles_ num_quantiles,
    std::vector<std::vector<Rcpp::NumericVector> >& quantiles,
    const bool compute_min_rank,
    std::vector<Rcpp::IntegerVector>& min_rank
) {
    ptrs.resize(num_groups);

    if (compute_min) {
        min.reserve(num_groups);
    }
    if (compute_mean) {
        mean.reserve(num_groups);
    }
    if (compute_median) {
        median.reserve(num_groups);
    }
    if (compute_max) {
        max.reserve(num_groups);
    }
    if (num_quantiles) {
        quantiles.reserve(num_groups);
    }
    if (compute_min_rank) {
        min_rank.reserve(num_groups);
    }

    for (I<decltype(num_groups)> g = 0; g < num_groups; ++g) {
        auto& curptr = ptrs[g];
        if (compute_min) {
            min.emplace_back(num_genes);
            curptr.min = min.back().begin();
        }
        if (compute_mean) {
            mean.emplace_back(num_genes);
            curptr.mean = mean.back().begin();
        }
        if (compute_median) {
            median.emplace_back(num_genes);
            curptr.median = median.back().begin();
        }
        if (compute_max) {
            max.emplace_back(num_genes);
            curptr.max = max.back().begin();
        }
        if (num_quantiles) {
            quantiles.emplace_back();
            quantiles.back().reserve(num_quantiles);
            curptr.quantiles.emplace();
            curptr.quantiles->reserve(num_quantiles);
            for (I<decltype(num_groups)> q = 0; q < num_quantiles; ++q) {
                quantiles.back().emplace_back(num_genes);
                curptr.quantiles->push_back(quantiles.back().back().begin());
            }
        }
        if (compute_min_rank) {
            min_rank.emplace_back(num_genes);
            curptr.min_rank = min_rank.back().begin();
        }
    }
}

inline auto setup_quantile_options(const Rcpp::Nullable<Rcpp::NumericVector>& input, std::optional<std::vector<double> >& output) {
    if (input.isNull()) {
        return sanisizer::as_size_type<Rcpp::NumericVector>(0);
    }
    Rcpp::NumericVector squantiles(input);
    output.emplace(squantiles.begin(), squantiles.end());
    return squantiles.size();
}

template<typename Ngroups_>
Rcpp::List format_summary_output(
    const Ngroups_ num_groups,
    const bool compute_min,
    std::vector<Rcpp::NumericVector>& min,
    const bool compute_mean,
    std::vector<Rcpp::NumericVector>& mean,
    const bool compute_median,
    std::vector<Rcpp::NumericVector>& median,
    const bool compute_max,
    std::vector<Rcpp::NumericVector>& max,
    const bool compute_quantiles,
    std::vector<std::vector<Rcpp::NumericVector> >& output_quantiles,
    const bool compute_min_rank,
    std::vector<Rcpp::IntegerVector>& min_rank
) { 
    auto output = sanisizer::create<Rcpp::List>(num_groups);
    for (I<decltype(num_groups)> g = 0; g < num_groups; ++g) {
        Rcpp::List current;

        if (compute_min) {
            current["min"] = std::move(min[g]);
        }
        if (compute_mean) {
            current["mean"] = std::move(mean[g]);
        }
        if (compute_median) {
            current["median"] = std::move(median[g]);
        }
        if (compute_max) {
            current["max"] = std::move(max[g]);
        }

        if (compute_quantiles) {
            auto& oquantiles = output_quantiles[g];
            const auto num_quantiles = oquantiles.size();
            auto collected = sanisizer::create<Rcpp::List>(num_quantiles);
            for (I<decltype(num_quantiles)> q = 0; q < num_quantiles; ++q) {
                collected[q] = std::move(oquantiles[q]);
            }
            current["quantile"] = std::move(collected);
        }

        if (compute_min_rank) {
            current["min.rank"] = std::move(min_rank[g]);
        }

        output[g] = std::move(current);
    }

    return output;
}

#endif
