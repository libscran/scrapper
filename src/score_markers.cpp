//#include "config.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "Rcpp.h"
#include "Rtatami.h"

#include "scran_markers/scran_markers.hpp"
#include "scran_markers/score_markers_best.hpp"

#include "utils_block.h"

void configure_group_vectors(Rcpp::NumericMatrix& store, std::vector<double*>& ptrs, int NR, int num_groups) { 
    store = Rcpp::NumericMatrix(NR, num_groups);
    ptrs.reserve(num_groups);
    for (int g = 0; g < num_groups; ++g) {
        const auto out_offset = sanisizer::product_unsafe<std::size_t>(g, NR);
        ptrs.emplace_back(store.begin() + out_offset);
    }
}

scran_markers::BlockAveragePolicy process_average_policy(const std::string& block_average_policy) {
    if (block_average_policy == "mean") {
        return scran_markers::BlockAveragePolicy::MEAN;
    } else if (block_average_policy == "quantile") {
        return scran_markers::BlockAveragePolicy::QUANTILE;
    } else {
        throw std::runtime_error("block average policy should be either 'mean' or 'quantile'");
        return scran_markers::BlockAveragePolicy::MEAN;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_summary(
    SEXP x,
    Rcpp::IntegerVector groups,
    int num_groups,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    std::string block_average_policy,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double block_quantile,
    double threshold,
    int num_threads,
    bool compute_group_mean,
    bool compute_group_detected,
    bool compute_delta_mean,
    bool compute_delta_detected,
    bool compute_cohens_d,
    bool compute_auc,
    int min_rank_limit
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersSummaryOptions opt;
    opt.threshold = threshold;
    opt.num_threads = num_threads;
    opt.min_rank_limit = min_rank_limit;
    opt.block_average_policy = process_average_policy(block_average_policy);
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.block_quantile = block_quantile;

    scran_markers::ScoreMarkersSummaryBuffers<double, int> buffers;

    Rcpp::NumericMatrix means, detected;
    if (compute_group_mean) {
        configure_group_vectors(means, buffers.mean, NR, num_groups);
    }
    if (compute_group_detected) {
        configure_group_vectors(detected, buffers.detected, NR, num_groups);
    }

    std::vector<Rcpp::NumericVector> cohens_min, cohens_mean, cohens_median, cohens_max;
    std::vector<Rcpp::NumericVector> auc_min, auc_mean, auc_median, auc_max;
    std::vector<Rcpp::NumericVector> dm_min, dm_mean, dm_median, dm_max;
    std::vector<Rcpp::NumericVector> dd_min, dd_mean, dd_median, dd_max;
    std::vector<Rcpp::IntegerVector> cohens_mr, auc_mr, dm_mr, dd_mr;

    const auto initialize = [&](
        std::vector<scran_markers::SummaryBuffers<double, int> >& ptrs,
        std::vector<Rcpp::NumericVector>& min,
        std::vector<Rcpp::NumericVector>& mean,
        std::vector<Rcpp::NumericVector>& median,
        std::vector<Rcpp::NumericVector>& max,
        std::vector<Rcpp::IntegerVector>& min_rank
    ) -> void {
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
        if (!sanisizer::is_equal(block_info.size(), NC)) {
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
    std::string block_average_policy,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double block_quantile,
    double threshold,
    int num_threads,
    bool compute_group_mean,
    bool compute_group_detected,
    bool compute_delta_mean,
    bool compute_delta_detected,
    bool compute_cohens_d,
    bool compute_auc)
{
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersPairwiseOptions opt;
    opt.threshold = threshold;
    opt.num_threads = num_threads;
    opt.block_average_policy = process_average_policy(block_average_policy);
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.block_quantile = block_quantile;

    scran_markers::ScoreMarkersPairwiseBuffers<double> buffers;

    Rcpp::NumericMatrix means, detected;
    if (compute_group_mean) {
        configure_group_vectors(means, buffers.mean, NR, num_groups);
    }
    if (compute_group_detected) {
        configure_group_vectors(detected, buffers.detected, NR, num_groups);
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
        if (!sanisizer::is_equal(block_info.size(), NC)) {
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

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_best(
    SEXP x,
    int top,
    Rcpp::IntegerVector groups,
    int num_groups,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    std::string block_average_policy,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double block_quantile,
    double threshold,
    int num_threads,
    bool compute_group_mean,
    bool compute_group_detected,
    bool compute_delta_mean,
    bool compute_delta_detected,
    bool compute_cohens_d,
    bool compute_auc)
{
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersBestOptions opt;
    opt.threshold = threshold;
    opt.num_threads = num_threads;
    opt.block_average_policy = process_average_policy(block_average_policy);
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.block_quantile = block_quantile;

    opt.compute_group_mean = compute_group_mean;
    opt.compute_group_detected = compute_group_detected;
    opt.compute_cohens_d = compute_cohens_d;
    opt.compute_auc = compute_auc;
    opt.compute_delta_mean = compute_delta_mean;
    opt.compute_delta_detected = compute_delta_detected;

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    scran_markers::ScoreMarkersBestResults<double, int> res; 
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), NC)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        res = scran_markers::score_markers_best_blocked<double>(*mat, static_cast<const int*>(groups.begin()), ptr, top, opt);
    } else {
        res = scran_markers::score_markers_best<double>(*mat, static_cast<const int*>(groups.begin()), top, opt);
    }

    const auto transfer_groupwise = [&](Rcpp::NumericMatrix& store, std::vector<std::vector<double> >& vecs) -> void {
        store = Rcpp::NumericMatrix(NR, num_groups);
        for (int g = 0; g < num_groups; ++g) {
            std::copy(vecs[g].begin(), vecs[g].end(), store.begin() + sanisizer::product_unsafe<std::size_t>(g, NR));
        }
    };
    Rcpp::NumericMatrix means, detected;
    if (compute_group_mean) {
        transfer_groupwise(means, res.mean);
    }
    if (compute_group_detected) {
        transfer_groupwise(detected, res.detected);
    }

    const auto transfer_effects = [&](Rcpp::List& store, std::vector<std::vector<std::vector<std::pair<int, double> > > >& vecs) -> void {
        store = Rcpp::List(num_groups);
        for (int g = 0; g < num_groups; ++g) {
            Rcpp::List current(num_groups);
            for (int g2 = 0; g2 < num_groups; ++g2) {
                if (g == g2) {
                    continue;
                }
                const auto& curtop = vecs[g][g2];
                const std::size_t numtop = curtop.size();
                Rcpp::IntegerVector indices(numtop);
                Rcpp::NumericVector effects(numtop);
                for (std::size_t t = 0; t < numtop; ++t) {
                    indices[t] = curtop[t].first + 1;
                    effects[t] = curtop[t].second;
                }
                current[g2] = Rcpp::DataFrame::create(
                    Rcpp::Named("index") = indices,
                    Rcpp::Named("effect") = effects
                );
            }
            store[g] = current;
        }
    };
    Rcpp::List cohens_d, auc, delta_mean, delta_detected;
    if (compute_cohens_d) {
        transfer_effects(cohens_d, res.cohens_d);
    }
    if (compute_auc) {
        transfer_effects(auc, res.auc);
    }
    if (compute_delta_mean) {
        transfer_effects(delta_mean, res.delta_mean);
    }
    if (compute_delta_detected) {
        transfer_effects(delta_detected, res.delta_detected);
    }

    return Rcpp::List::create(
        Rcpp::Named("mean") = means,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("cohens.d") = cohens_d,
        Rcpp::Named("auc") = auc,
        Rcpp::Named("delta.mean") = delta_mean,
        Rcpp::Named("delta.detected") = delta_detected
    );
}
