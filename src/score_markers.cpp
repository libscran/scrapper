//#include "config.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "Rcpp.h"
#include "Rtatami.h"

#include "scran_markers/scran_markers.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_summary(
    SEXP x,
    Rcpp::IntegerVector groups,
    int num_groups,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double threshold,
    int num_threads,
    bool compute_delta_mean,
    bool compute_delta_detected,
    bool compute_cohens_d,
    bool compute_auc)
{
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    size_t NC = mat->ncol();
    size_t NR = mat->nrow();
    if (static_cast<size_t>(groups.size()) != NC) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.threshold = threshold;
    opt.num_threads = num_threads;
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);

    scran_markers::ScoreMarkersSummaryBuffers<double, int> buffers;

    Rcpp::NumericMatrix means(NR, num_groups);
    Rcpp::NumericMatrix detected(NR, num_groups);
    buffers.mean.reserve(num_groups);
    buffers.detected.reserve(num_groups);
    size_t out_offset = 0;
    for (int g = 0; g < num_groups; ++g) {
        buffers.mean.emplace_back(means.begin() + out_offset);
        buffers.detected.emplace_back(detected.begin() + out_offset);
        out_offset += NR;
    }

    std::vector<Rcpp::NumericVector> cohens_min, cohens_mean, cohens_median, cohens_max;
    std::vector<Rcpp::NumericVector> auc_min, auc_mean, auc_median, auc_max;
    std::vector<Rcpp::NumericVector> dm_min, dm_mean, dm_median, dm_max;
    std::vector<Rcpp::NumericVector> dd_min, dd_mean, dd_median, dd_max;
    std::vector<Rcpp::IntegerVector> cohens_mr, auc_mr, dm_mr, dd_mr;

    auto initialize = [&](
        std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs,
        std::vector<Rcpp::NumericVector>& min,
        std::vector<Rcpp::NumericVector>& mean,
        std::vector<Rcpp::NumericVector>& median,
        std::vector<Rcpp::NumericVector>& max,
        std::vector<Rcpp::IntegerVector>& min_rank
    ) {
        ptrs.resize(num_groups);
        min.reserve(num_groups);
        mean.reserve(num_groups);
        median.reserve(num_groups);
        max.reserve(num_groups);
        min_rank.reserve(num_groups);
        for (int g = 0; g < num_groups; ++g) {
            min.emplace_back(NR);
            ptrs[g].min = min.back().begin();
            mean.emplace_back(NR);
            ptrs[g].mean = mean.back().begin();
            median.emplace_back(NR);
            ptrs[g].median = median.back().begin();
            max.emplace_back(NR);
            ptrs[g].max = max.back().begin();
            min_rank.emplace_back(NR);
            ptrs[g].min_rank = min_rank.back().begin();
        }
    };

    if (compute_cohens_d) {
        initialize(buffers.cohens_d, cohens_min, cohens_mean, cohens_median, cohens_max, cohens_mr);
    }
    if (compute_delta_mean) {
        initialize(buffers.delta_mean, dm_min, dm_mean, dm_median, dm_max, dm_mr);
    }
    if (compute_delta_detected) {
        initialize(buffers.delta_detected, dd_min, dd_mean, dd_median, dd_max, dd_mr);
    }
    if (compute_auc) {
        initialize(buffers.auc, auc_min, auc_mean, auc_median, auc_max, auc_mr);
    }

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (block_info.size() != NC) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        scran_markers::score_markers_summary_blocked(*mat, static_cast<const int*>(groups.begin()), ptr, opt, buffers);
    } else {
        scran_markers::score_markers_summary(*mat, static_cast<const int*>(groups.begin()), opt, buffers);
    }

    auto format = [&](
        const std::vector<Rcpp::NumericVector>& min,
        const std::vector<Rcpp::NumericVector>& mean,
        const std::vector<Rcpp::NumericVector>& median,
        const std::vector<Rcpp::NumericVector>& max,
        const std::vector<Rcpp::IntegerVector>& min_rank
    ) -> Rcpp::List {
        size_t ngroups = min.size();
        Rcpp::List output(ngroups);
        for (size_t g = 0; g < ngroups; ++g) {
            output[g] = Rcpp::List::create(
                Rcpp::Named("min") = min[g],
                Rcpp::Named("mean") = mean[g],
                Rcpp::Named("median") = median[g],
                Rcpp::Named("max") = max[g],
                Rcpp::Named("min.rank") = min_rank[g]
            );
        }
        return output;
    };

    return Rcpp::List::create(
        Rcpp::Named("mean") = means,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("cohens.d") = format(cohens_min, cohens_mean, cohens_median, cohens_max, cohens_mr),
        Rcpp::Named("auc") = format(auc_min, auc_mean, auc_median, auc_max, auc_mr),
        Rcpp::Named("delta.mean") = format(dm_min, dm_mean, dm_median, dm_max, dm_mr),
        Rcpp::Named("delta.detected") = format(dd_min, dd_mean, dd_median, dd_max, dd_mr)
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_pairwise(
    SEXP x,
    Rcpp::IntegerVector groups,
    int num_groups,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double threshold,
    int num_threads,
    bool compute_delta_mean,
    bool compute_delta_detected,
    bool compute_cohens_d,
    bool compute_auc)
{
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    size_t NC = mat->ncol();
    size_t NR = mat->nrow();
    if (static_cast<size_t>(groups.size()) != NC) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.threshold = threshold;
    opt.num_threads = num_threads;
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);

    scran_markers::ScoreMarkersPairwiseBuffers<double> buffers;

    Rcpp::NumericMatrix means(NR, num_groups);
    Rcpp::NumericMatrix detected(NR, num_groups);
    buffers.mean.reserve(num_groups);
    buffers.detected.reserve(num_groups);
    size_t out_offset = 0;
    for (int g = 0; g < num_groups; ++g) {
        buffers.mean.emplace_back(means.begin() + out_offset);
        buffers.detected.emplace_back(detected.begin() + out_offset);
        out_offset += NR;
    }

    Rcpp::Dimension dim(num_groups, num_groups, NR);
    Rcpp::NumericVector cohen, auc, delta_mean, delta_detected;
    if (compute_cohens_d) {
        cohen = Rcpp::NumericVector(dim);
        buffers.cohens_d = cohen.begin();
    }
    if (compute_delta_mean) {
        delta_mean = Rcpp::NumericVector(dim);
        buffers.delta_mean = delta_mean.begin();
    }
    if (compute_delta_detected) {
        delta_detected = Rcpp::NumericVector(dim);
        buffers.delta_detected = delta_detected.begin();
    }
    if (compute_auc) {
        auc = Rcpp::NumericVector(dim);
        buffers.auc = auc.begin();
    }

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (block_info.size() != NC) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        scran_markers::score_markers_pairwise_blocked(*mat, static_cast<const int*>(groups.begin()), ptr, opt, buffers);
    } else {
        scran_markers::score_markers_pairwise(*mat, static_cast<const int*>(groups.begin()), opt, buffers);
    }

    return Rcpp::List::create(
        Rcpp::Named("mean") = means,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("cohens.d") = cohen,
        Rcpp::Named("auc") = auc,
        Rcpp::Named("delta.mean") = delta_mean,
        Rcpp::Named("delta.detected") = delta_detected
    );
}
