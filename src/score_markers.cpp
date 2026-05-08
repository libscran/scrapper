#include "config.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "scran_markers/scran_markers.hpp"
#include "scran_markers/score_markers_best.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_markers.h"
#include "utils_other.h"

template<typename Nrow_, typename Ngroups_>
void configure_group_vectors(Rcpp::NumericMatrix& store, std::vector<double*>& ptrs, Nrow_ NR, Ngroups_ num_groups) { 
    store = create_matrix<Rcpp::NumericMatrix>(NR, num_groups);
    sanisizer::reserve(ptrs, num_groups);
    for (I<decltype(num_groups)> g = 0; g < num_groups; ++g) {
        const auto out_offset = sanisizer::product_unsafe<std::size_t>(g, NR);
        ptrs.emplace_back(store.begin() + out_offset);
    }
}

void set_block_average_policy(const Rcpp::RObject& block_average_policy, scran_markers::BlockAveragePolicy& target) {
    if (block_average_policy.isNULL()) {
        return;
    }
    const auto bap = parse_single_string(block_average_policy, "block.average.policy");
    if (bap == "mean") {
        target = scran_markers::BlockAveragePolicy::MEAN;
    } else if (bap == "quantile") {
        target = scran_markers::BlockAveragePolicy::QUANTILE;
    } else {
        throw std::runtime_error("block.average.policy should be either 'mean' or 'quantile'");
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_summary(
    SEXP x,
    Rcpp::IntegerVector groups,
    int num_groups,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    Rcpp::RObject block_average_policy,
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject block_quantile,
    Rcpp::RObject threshold,
    Rcpp::RObject num_threads,
    Rcpp::RObject compute_group_mean,
    Rcpp::RObject compute_group_detected,
    Rcpp::RObject compute_delta_mean,
    Rcpp::RObject compute_delta_detected,
    Rcpp::RObject compute_cohens_d,
    Rcpp::RObject compute_auc,
    Rcpp::RObject compute_summary_min,
    Rcpp::RObject compute_summary_mean,
    Rcpp::RObject compute_summary_median,
    Rcpp::RObject compute_summary_max,
    Rcpp::RObject compute_summary_quantiles,
    Rcpp::RObject compute_summary_min_rank,
    Rcpp::RObject min_rank_limit
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersSummaryOptions opt;
    set_block_average_policy(block_average_policy, opt.block_average_policy);
    set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
    set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");
    set_number(block_quantile, opt.block_quantile, "block.quantile");
    set_number(threshold, opt.threshold, "threshold");
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_bool(compute_group_mean, opt.compute_group_mean, "compute.group.mean"); // technically these are not necessary, but we'll set it for consistency.
    set_bool(compute_group_detected, opt.compute_group_detected, "compute.group.detected");
    set_bool(compute_delta_mean, opt.compute_delta_mean, "compute.delta.mean");
    set_bool(compute_delta_detected, opt.compute_delta_detected, "compute.delta.detected");
    set_bool(compute_cohens_d, opt.compute_cohens_d, "compute.cohens.d");
    set_bool(compute_auc, opt.compute_auc, "compute.auc");
    set_bool(compute_summary_min, opt.compute_min, "compute.summary.min");
    set_bool(compute_summary_mean, opt.compute_mean, "compute.summary.mean");
    set_bool(compute_summary_median, opt.compute_median, "compute.summary.median");
    set_bool(compute_summary_max, opt.compute_max, "compute.summary.max");
    const auto num_quantiles = setup_summary_quantiles(compute_summary_quantiles, opt.compute_summary_quantiles);
    set_bool(compute_summary_min_rank, opt.compute_min_rank, "compute.summary.min.rank");
    set_integer(min_rank_limit, opt.min_rank_limit, "min.rank.limit");

    scran_markers::ScoreMarkersSummaryBuffers<double, int> buffers;

    Rcpp::NumericMatrix means, detected;
    if (opt.compute_group_mean) {
        configure_group_vectors(means, buffers.mean, NR, num_groups);
    }
    if (opt.compute_group_detected) {
        configure_group_vectors(detected, buffers.detected, NR, num_groups);
    }

    std::vector<Rcpp::NumericVector> cohens_min, cohens_mean, cohens_median, cohens_max;
    std::vector<Rcpp::NumericVector> auc_min, auc_mean, auc_median, auc_max;
    std::vector<Rcpp::NumericVector> dm_min, dm_mean, dm_median, dm_max;
    std::vector<Rcpp::NumericVector> dd_min, dd_mean, dd_median, dd_max;
    std::vector<std::vector<Rcpp::NumericVector> > cohens_quant, auc_quant, dm_quant, dd_quant;
    std::vector<Rcpp::IntegerVector> cohens_mr, auc_mr, dm_mr, dd_mr;

    if (opt.compute_cohens_d) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.cohens_d,
            opt.compute_min,
            cohens_min,
            opt.compute_mean,
            cohens_mean,
            opt.compute_median,
            cohens_median,
            opt.compute_max,
            cohens_max,
            num_quantiles,
            cohens_quant,
            opt.compute_min_rank,
            cohens_mr
        );
    }

    if (opt.compute_auc) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.auc,
            opt.compute_min,
            auc_min,
            opt.compute_mean,
            auc_mean,
            opt.compute_median,
            auc_median,
            opt.compute_max,
            auc_max,
            num_quantiles,
            auc_quant,
            opt.compute_min_rank,
            auc_mr
        );
    }

    if (opt.compute_delta_mean) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.delta_mean,
            opt.compute_min,
            dm_min,
            opt.compute_mean,
            dm_mean,
            opt.compute_median,
            dm_median,
            opt.compute_max,
            dm_max,
            num_quantiles,
            dm_quant,
            opt.compute_min_rank,
            dm_mr
        );
    }

    if (opt.compute_delta_detected) {
        initialize_summary_buffers(
            num_groups,
            NR,
            buffers.delta_detected,
            opt.compute_min,
            dd_min,
            opt.compute_mean,
            dd_mean,
            opt.compute_median,
            dd_median,
            opt.compute_max,
            dd_max,
            num_quantiles,
            dd_quant,
            opt.compute_min_rank,
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
    if (opt.compute_group_mean) {
        output["mean"] = std::move(means);
    }
    if (opt.compute_group_detected) {
        output["detected"] = std::move(detected);
    }

    if (opt.compute_cohens_d) {
        output["cohens.d"] = format_summary_output(
            num_groups,
            opt.compute_min,
            cohens_min,
            opt.compute_mean,
            cohens_mean,
            opt.compute_median,
            cohens_median,
            opt.compute_max,
            cohens_max,
            opt.compute_summary_quantiles.has_value(),
            cohens_quant,
            opt.compute_min_rank,
            cohens_mr
        );
    }

    if (opt.compute_auc) {
        output["auc"] = format_summary_output(
            num_groups,
            opt.compute_min,
            auc_min,
            opt.compute_mean,
            auc_mean,
            opt.compute_median,
            auc_median,
            opt.compute_max,
            auc_max,
            opt.compute_summary_quantiles.has_value(),
            auc_quant, 
            opt.compute_min_rank,
            auc_mr
        );
    }

    if (opt.compute_delta_mean) {
        output["delta.mean"] = format_summary_output(
            num_groups,
            opt.compute_min,
            dm_min,
            opt.compute_mean,
            dm_mean,
            opt.compute_median,
            dm_median,
            opt.compute_max,
            dm_max,
            opt.compute_summary_quantiles.has_value(),
            dm_quant,
            opt.compute_min_rank,
            dm_mr
        );
    }

    if (opt.compute_delta_detected) {
        output["delta.detected"] = format_summary_output(
            num_groups,
            opt.compute_min,
            dd_min,
            opt.compute_mean,
            dd_mean,
            opt.compute_median,
            dd_median,
            opt.compute_max,
            dd_max,
            opt.compute_summary_quantiles.has_value(),
            dd_quant,
            opt.compute_min_rank,
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
    Rcpp::RObject block_average_policy,
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject block_quantile,
    Rcpp::RObject threshold,
    Rcpp::RObject num_threads,
    Rcpp::RObject compute_group_mean,
    Rcpp::RObject compute_group_detected,
    Rcpp::RObject compute_delta_mean,
    Rcpp::RObject compute_delta_detected,
    Rcpp::RObject compute_cohens_d,
    Rcpp::RObject compute_auc
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersPairwiseOptions opt;
    set_block_average_policy(block_average_policy, opt.block_average_policy);
    set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
    set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");
    set_number(block_quantile, opt.block_quantile, "block.quantile");
    set_number(threshold, opt.threshold, "threshold");
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_bool(compute_group_mean, opt.compute_group_mean, "compute.group.mean"); // technically these are not necessary, but we'll set it for consistency.
    set_bool(compute_group_detected, opt.compute_group_detected, "compute.group.detected");
    set_bool(compute_delta_mean, opt.compute_delta_mean, "compute.delta.mean");
    set_bool(compute_delta_detected, opt.compute_delta_detected, "compute.delta.detected");
    set_bool(compute_cohens_d, opt.compute_cohens_d, "compute.cohens.d");
    set_bool(compute_auc, opt.compute_auc, "compute.auc");

    scran_markers::ScoreMarkersPairwiseBuffers<double> buffers;

    Rcpp::NumericMatrix means, detected;
    if (opt.compute_group_mean) {
        configure_group_vectors(means, buffers.mean, NR, num_groups);
    }
    if (opt.compute_group_detected) {
        configure_group_vectors(detected, buffers.detected, NR, num_groups);
    }

    Rcpp::Dimension dim(
        sanisizer::cast<std::size_t>(num_groups),
        num_groups, // also a size_t, but safety is covered by the previous statement.
        sanisizer::cast<std::size_t>(NR)
    );

    Rcpp::NumericVector cohens_d, auc, delta_mean, delta_detected;
    if (opt.compute_cohens_d) {
        cohens_d = Rcpp::NumericVector(dim);
        buffers.cohens_d = cohens_d.begin();
    }
    if (opt.compute_delta_mean) {
        delta_mean = Rcpp::NumericVector(dim);
        buffers.delta_mean = delta_mean.begin();
    }
    if (opt.compute_delta_detected) {
        delta_detected = Rcpp::NumericVector(dim);
        buffers.delta_detected = delta_detected.begin();
    }
    if (opt.compute_auc) {
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
    if (opt.compute_group_mean) {
        output["mean"] = std::move(means);
    }
    if (opt.compute_group_detected) {
        output["detected"] = std::move(detected);
    }

    if (opt.compute_cohens_d) {
        output["cohens.d"] = std::move(cohens_d);
    }
    if (opt.compute_auc) {
        output["auc"] = std::move(auc);
    }
    if (opt.compute_delta_mean) {
        output["delta.mean"] = std::move(delta_mean);
    }
    if (opt.compute_delta_detected) {
        output["delta.detected"] = std::move(delta_detected);
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
    Rcpp::RObject block_average_policy,
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject block_quantile,
    Rcpp::RObject threshold,
    Rcpp::RObject num_threads,
    Rcpp::RObject compute_group_mean,
    Rcpp::RObject compute_group_detected,
    Rcpp::RObject compute_delta_mean,
    Rcpp::RObject compute_delta_detected,
    Rcpp::RObject compute_cohens_d,
    Rcpp::RObject compute_auc,
    bool index_only
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();
    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("'groups' must have length equal to the number of cells");
    }

    scran_markers::ScoreMarkersBestOptions opt;
    set_block_average_policy(block_average_policy, opt.block_average_policy);
    set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
    set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");
    set_number(threshold, opt.threshold, "threshold");
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_number(block_quantile, opt.block_quantile, "block.quantile");
    set_bool(compute_group_mean, opt.compute_group_mean, "compute.group.mean");
    set_bool(compute_group_detected, opt.compute_group_detected, "compute.group.detected");
    set_bool(compute_delta_mean, opt.compute_delta_mean, "compute.delta.mean");
    set_bool(compute_delta_detected, opt.compute_delta_detected, "compute.delta.detected");
    set_bool(compute_cohens_d, opt.compute_cohens_d, "compute.cohens.d");
    set_bool(compute_auc, opt.compute_auc, "compute.auc");

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
        for (I<decltype(num_groups)> g = 0; g < num_groups; ++g) {
            std::copy(vecs[g].begin(), vecs[g].end(), store.begin() + sanisizer::product_unsafe<std::size_t>(g, NR));
        }
    };
    Rcpp::NumericMatrix means, detected;
    if (opt.compute_group_mean) {
        transfer_groupwise(means, res.mean);
    }
    if (opt.compute_group_detected) {
        transfer_groupwise(detected, res.detected);
    }

    const auto transfer_effects = [&](Rcpp::List& store, std::vector<std::vector<std::vector<std::pair<int, double> > > >& vecs) -> void {
        store = sanisizer::create<Rcpp::List>(num_groups);
        for (I<decltype(num_groups)> g = 0; g < num_groups; ++g) {
            auto current = sanisizer::create<Rcpp::List>(num_groups);
            for (I<decltype(num_groups)> g2 = 0; g2 < num_groups; ++g2) {
                if (g == g2) {
                    continue;
                }

                const auto& curtop = vecs[g][g2];
                const auto numtop = curtop.size();
                auto indices = sanisizer::create<Rcpp::IntegerVector>(numtop);

                if (index_only) {
                    for (I<decltype(numtop)> t = 0; t < numtop; ++t) {
                        indices[t] = curtop[t].first + 1;
                    }
                    current[g2] = std::move(indices);

                } else {
                    auto effects = sanisizer::create<Rcpp::NumericVector>(numtop);
                    for (I<decltype(numtop)> t = 0; t < numtop; ++t) {
                        indices[t] = curtop[t].first + 1;
                        effects[t] = curtop[t].second;
                    }

                    Rcpp::S4 df("DFrame");
                    df.slot("rownames") = R_NilValue;
                    df.slot("nrows") = sanisizer::cast<int>(numtop);
                    df.slot("listData") = Rcpp::List::create(
                        Rcpp::Named("index") = std::move(indices),
                        Rcpp::Named("effect") = std::move(effects)
                    );
                    df.slot("elementType") = "ANY";
                    df.slot("elementMetadata") = R_NilValue;
                    df.slot("metadata") = Rcpp::List();
                    current[g2] = std::move(df);
                }
            }
            store[g] = std::move(current);
        }
    };

    Rcpp::List cohens_d, auc, delta_mean, delta_detected;
    if (opt.compute_cohens_d) {
        transfer_effects(cohens_d, res.cohens_d);
    }
    if (opt.compute_auc) {
        transfer_effects(auc, res.auc);
    }
    if (opt.compute_delta_mean) {
        transfer_effects(delta_mean, res.delta_mean);
    }
    if (opt.compute_delta_detected) {
        transfer_effects(delta_detected, res.delta_detected);
    }

    Rcpp::List output;
    if (opt.compute_group_mean) {
        output["mean"] = std::move(means);
    }
    if (opt.compute_group_detected) {
        output["detected"] = std::move(detected);
    }

    if (opt.compute_cohens_d) {
        output["cohens.d"] = std::move(cohens_d);
    }
    if (opt.compute_auc) {
        output["auc"] = std::move(auc);
    }
    if (opt.compute_delta_mean) {
        output["delta.mean"] = std::move(delta_mean);
    }
    if (opt.compute_delta_detected) {
        output["delta.detected"] = std::move(delta_detected);
    }

    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List score_markers_defaults(int mode) {
    Rcpp::List output;

    auto format_common_defaults = [&](const auto& opt) -> void {
        if (opt.block_average_policy == scran_markers::BlockAveragePolicy::MEAN) {
            output["block.average.policy"] = "mean";
        } else {
            throw std::runtime_error("unexpected block.average.policy default for scoreMarkers");
        }

        report_block_weight_policy_default(output, opt.block_weight_policy, "block.weight.policy", "scoreMarkers");
        report_variable_block_weight_default(output, opt.variable_block_weight_parameters, "variable.block.weight");
        output["block.quantile"] = opt.block_quantile;
        output["threshold"] = opt.threshold;
        output["num.threads"] = opt.num_threads;
        output["compute.group.mean"] = opt.compute_group_mean;
        output["compute.group.detected"] = opt.compute_group_detected;
        output["compute.delta.mean"] = opt.compute_delta_mean;
        output["compute.delta.detected"] = opt.compute_delta_detected;
        output["compute.cohens.d"] = opt.compute_cohens_d;
        output["compute.auc"] = opt.compute_auc;
    };

    if (mode == 0) {
        scran_markers::ScoreMarkersSummaryOptions opt;
        format_common_defaults(opt);
    } else if (mode == 1) {
        scran_markers::ScoreMarkersPairwiseOptions opt;
        format_common_defaults(opt);
    } else {
        scran_markers::ScoreMarkersBestOptions opt;
        format_common_defaults(opt);
    }

    // Adding all the summary options.
    scran_markers::ScoreMarkersSummaryOptions opt;
    output["compute.summary.min"] = opt.compute_min;
    output["compute.summary.mean"] = opt.compute_mean;
    output["compute.summary.median"] = opt.compute_median;
    output["compute.summary.max"] = opt.compute_max;
    if (opt.compute_summary_quantiles.has_value()) {
        throw std::runtime_error("unexpected compute.summary.quantiles default for scoreMarkers");
    } else {
        output["compute.summary.quantiles"] = R_NilValue;
    }
    output["compute.summary.min.rank"] = opt.compute_min_rank;
    output["min.rank.limit"] = opt.min_rank_limit;

    return output;
}
