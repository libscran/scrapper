#ifndef UTILS_MARKERS_H
#define UTILS_MARKERS_H

#include "Rcpp.h"

#include <vector>
#include <cstddef>
#include <string>

#include "scran_markers/scran_markers.hpp"

inline void initialize_summary_buffers(
    const int num_groups,
    const int num_genes,
    std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs,
    const bool compute_min,
    std::vector<Rcpp::NumericVector>& min,
    const bool compute_mean,
    std::vector<Rcpp::NumericVector>& mean,
    const bool compute_median,
    std::vector<Rcpp::NumericVector>& median,
    const bool compute_max,
    std::vector<Rcpp::NumericVector>& max,
    const std::size_t num_quantiles,
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

    for (int g = 0; g < num_groups; ++g) {
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
            for (std::size_t q = 0; q < num_quantiles; ++q) {
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

inline std::size_t setup_quantile_options(const Rcpp::Nullable<Rcpp::NumericVector>& input, std::optional<std::vector<double> >& output) {
    if (input.isNull()) {
        return 0;
    }
    Rcpp::NumericVector squantiles(input);
    output.emplace(squantiles.begin(), squantiles.end());
    return sanisizer::cast<std::size_t>(squantiles.size());
}

inline Rcpp::List format_summary_output(
    const int num_groups,
    const bool compute_min,
    const std::vector<Rcpp::NumericVector>& min,
    const bool compute_mean,
    const std::vector<Rcpp::NumericVector>& mean,
    const bool compute_median,
    const std::vector<Rcpp::NumericVector>& median,
    const bool compute_max,
    const std::vector<Rcpp::NumericVector>& max,
    const bool compute_quantiles,
    const std::vector<std::vector<Rcpp::NumericVector> >& output_quantiles,
    const bool compute_min_rank,
    const std::vector<Rcpp::IntegerVector>& min_rank
) { 
    auto output = sanisizer::create<Rcpp::List>(num_groups);
    for (int g = 0; g < num_groups; ++g) {
        Rcpp::List current;
        if (compute_min) {
            current["min"] = min[g];
        }
        if (compute_mean) {
            current["mean"] = mean[g];
        }
        if (compute_median) {
            current["median"] = median[g];
        }
        if (compute_max) {
            current["max"] = max[g];
        }
        if (compute_quantiles) {
            const auto& oquantiles = output_quantiles[g];
            auto num_quantiles = oquantiles.size();
            auto collected = sanisizer::create<Rcpp::List>(num_quantiles);
            for (decltype(num_quantiles) q = 0; q < num_quantiles; ++q) {
                collected[q] = oquantiles[q];
            }
            current["quantile"] = std::move(collected);
        }
        if (compute_min_rank) {
            current["min.rank"] = min_rank[g];
        }
        output[g] = current;
    }
    return output;
}

#endif
