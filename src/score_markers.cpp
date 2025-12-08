//#include "config.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "Rcpp.h"
#include "Rtatami.h"

#include "scran_markers/scran_markers.hpp"
#include "scran_markers/score_markers_best.hpp"

#include "utils_block.h"
#include "utils_markers.h"

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
    bool compute_summary_min,
    bool compute_summary_mean,
    bool compute_summary_median,
    bool compute_summary_max,
    Rcpp::Nullable<Rcpp::NumericVector> compute_summary_quantiles,
    bool compute_summary_min_rank,
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
    const std::size_t num_quantiles = setup_quantile_options(compute_summary_quantiles, opt.compute_summary_quantiles);

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
    std::vector<std::vector<Rcpp::NumericVector> > cohens_quant, auc_quant, dm_quant, dd_quant;
    std::vector<Rcpp::IntegerVector> cohens_mr, auc_mr, dm_mr, dd_mr;

    if (compute_cohens_d) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.cohens_d,
            compute_summary_min,
            cohens_min,
            compute_summary_mean,
            cohens_mean,
            compute_summary_median,
            cohens_median,
            compute_summary_max,
            cohens_max,
            num_quantiles,
            cohens_quant,
            compute_summary_min_rank,
            cohens_mr
        );
    }

    if (compute_auc) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.auc,
            compute_summary_min,
            auc_min,
            compute_summary_mean,
            auc_mean,
            compute_summary_median,
            auc_median,
            compute_summary_max,
            auc_max,
            num_quantiles,
            auc_quant,
            compute_summary_min_rank,
            auc_mr
        );
    }

    if (compute_delta_mean) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.delta_mean,
            compute_summary_min,
            dm_min,
            compute_summary_mean,
            dm_mean,
            compute_summary_median,
            dm_median,
            compute_summary_max,
            dm_max,
            num_quantiles,
            dm_quant,
            compute_summary_min_rank,
            dm_mr
        );
    }

    if (compute_delta_detected) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.delta_detected,
            compute_summary_min,
            dd_min,
            compute_summary_mean,
            dd_mean,
            compute_summary_median,
            dd_median,
            compute_summary_max,
            dd_max,
            num_quantiles,
            dd_quant,
            compute_summary_min_rank,
            dd_mr
        );
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

    Rcpp::List output;
    if (compute_group_mean) {
        output["mean"] = means;
    }
    if (compute_group_detected) {
        output["detected"] = detected;
    }

    if (compute_cohens_d) {
        output["cohens.d"] = format_summary_output(
            num_groups,
            compute_summary_min,
            cohens_min,
            compute_summary_mean,
            cohens_mean,
            compute_summary_median,
            cohens_median,
            compute_summary_max,
            cohens_max,
            num_quantiles,
            opt.compute_summary_quantiles,
            cohens_quant,
            compute_summary_min_rank,
            cohens_mr
        );
    }

    if (compute_auc) {
        output["auc"] = format_summary_output(
            num_groups,
            compute_summary_min,
            auc_min,
            compute_summary_mean,
            auc_mean,
            compute_summary_median,
            auc_median,
            compute_summary_max,
            auc_max,
            num_quantiles,
            opt.compute_summary_quantiles,
            auc_quant, 
            compute_summary_min_rank,
            auc_mr
        );
    }

    if (compute_delta_mean) {
        output["delta.mean"] = format_summary_output(
            num_groups,
            compute_summary_min,
            dm_min,
            compute_summary_mean,
            dm_mean,
            compute_summary_median,
            dm_median,
            compute_summary_max,
            dm_max,
            num_quantiles,
            opt.compute_summary_quantiles,
            dm_quant,
            compute_summary_min_rank,
            dm_mr
        );
    }

    if (compute_delta_detected) {
        output["delta.detected"] = format_summary_output(
            num_groups,
            compute_summary_min,
            dd_min,
            compute_summary_mean,
            dd_mean,
            compute_summary_median,
            dd_median,
            compute_summary_max,
            dd_max,
            num_quantiles,
            opt.compute_summary_quantiles,
            dd_quant,
            compute_summary_min_rank,
            dd_mr
        );
    }

    return output;
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
    Rcpp::NumericVector cohens_d, auc, delta_mean, delta_detected;
    if (compute_cohens_d) {
        cohens_d = Rcpp::NumericVector(dim);
        buffers.cohens_d = cohens_d.begin();
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

    Rcpp::List output;
    if (compute_group_mean) {
        output["mean"] = means;
    }
    if (compute_group_detected) {
        output["detected"] = detected;
    }

    if (compute_cohens_d) {
        output["cohens.d"] = cohens_d;
    }
    if (compute_auc) {
        output["auc"] = auc;
    }
    if (compute_delta_mean) {
        output["delta.mean"] = delta_mean;
    }
    if (compute_delta_detected) {
        output["delta.detected"] = delta_detected;
    }

    return output;
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

    Rcpp::List output;
    if (compute_group_mean) {
        output["mean"] = means;
    }
    if (compute_group_detected) {
        output["detected"] = detected;
    }

    if (compute_cohens_d) {
        output["cohens.d"] = cohens_d;
    }
    if (compute_auc) {
        output["auc"] = auc;
    }
    if (compute_delta_mean) {
        output["delta.mean"] = delta_mean;
    }
    if (compute_delta_detected) {
        output["delta.detected"] = delta_detected;
    }

    return output;
}
