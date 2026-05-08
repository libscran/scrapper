#include "config.h"

#include <vector>
#include <stdexcept>

#include "scran_markers/scran_markers.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_markers.h"

//[[Rcpp::export(rng=false)]]
SEXP summarize_effects(
    int num_genes,
    int num_groups,
    Rcpp::NumericVector effects,
    Rcpp::RObject compute_summary_min,
    Rcpp::RObject compute_summary_mean,
    Rcpp::RObject compute_summary_median,
    Rcpp::RObject compute_summary_max,
    Rcpp::RObject compute_summary_quantiles,
    Rcpp::RObject compute_summary_min_rank,
    Rcpp::RObject num_threads
) {
    const auto expected = sanisizer::product<std::size_t>(num_groups, num_groups, num_genes);
    if (!sanisizer::is_equal(effects.size(), expected)) {
        throw std::runtime_error("'effects' does not have the expected length");
    }

    scran_markers::SummarizeEffectsOptions opt;
    set_bool(compute_summary_min, opt.compute_min, "compute.summary.min");
    set_bool(compute_summary_mean, opt.compute_mean, "compute.summary.mean");
    set_bool(compute_summary_median, opt.compute_median, "compute.summary.median");
    set_bool(compute_summary_median, opt.compute_median, "compute.summary.median");
    set_bool(compute_summary_max, opt.compute_max, "compute.summary.max");
    const auto num_quantiles = setup_summary_quantiles(compute_summary_quantiles, opt.compute_quantiles);
    set_bool(compute_summary_min_rank, opt.compute_min_rank, "compute.summary.min.rank");
    set_integer(num_threads, opt.num_threads, "num.threads");

    std::vector<Rcpp::NumericVector> min, mean, median, max;
    std::vector<std::vector<Rcpp::NumericVector> > quant;
    std::vector<Rcpp::IntegerVector> mr;

    std::vector<scran_markers::SummaryBuffers<double, int> > groupwise;
    initialize_summary_buffers(
        num_groups,
        num_genes,
        groupwise,
        opt.compute_min,
        min,
        opt.compute_mean,
        mean,
        opt.compute_median,
        median,
        opt.compute_max,
        max,
        num_quantiles,
        quant,
        opt.compute_min_rank,
        mr
    );

    scran_markers::summarize_effects(
        num_genes,
        num_groups,
        static_cast<const double*>(effects.begin()),
        groupwise,
        opt
    );

    return format_summary_output(
        num_groups,
        opt.compute_min,
        min,
        opt.compute_mean,
        mean,
        opt.compute_median,
        median,
        opt.compute_max,
        max,
        opt.compute_quantiles.has_value(),
        quant,
        opt.compute_min_rank,
        mr
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List summarize_effects_defaults() {
    Rcpp::List output;
    scran_markers::SummarizeEffectsOptions opt;
    output["compute.summary.min"] = opt.compute_min;
    output["compute.summary.mean"] = opt.compute_mean;
    output["compute.summary.median"] = opt.compute_median;
    output["compute.summary.max"] = opt.compute_max;
    if (opt.compute_quantiles.has_value()) {
        throw std::runtime_error("unexpected compute.summary.quantiles default for summarizeEffects");
    } else {
        output["compute.summary.quantiles"] = R_NilValue;
    }
    output["compute.summary.min.rank"] = opt.compute_min_rank;
    output["num.threads"] = opt.num_threads;
    return output;
}
