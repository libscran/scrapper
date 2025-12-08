//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"

#include "scran_markers/scran_markers.hpp"

#include "utils_markers.h"

//[[Rcpp::export(rng=false)]]
SEXP summarize_effects(
    int num_genes,
    int num_groups,
    Rcpp::NumericVector effects,
    bool compute_summary_min,
    bool compute_summary_mean,
    bool compute_summary_median,
    bool compute_summary_max,
    Rcpp::Nullable<Rcpp::NumericVector> compute_summary_quantiles,
    bool compute_summary_min_rank,
    int num_threads
) {
    const std::size_t expected = sanisizer::product<std::size_t>(num_groups, num_groups, num_genes);
    if (!sanisizer::is_equal(effects.size(), expected)) {
        throw std::runtime_error("'effects' does not have the expected length");
    }

    scran_markers::SummarizeEffectsOptions opt;
    opt.num_threads = num_threads;
    const std::size_t num_quantiles = setup_quantile_options(compute_summary_quantiles, opt.compute_quantiles);

    std::vector<Rcpp::NumericVector> min, mean, median, max;
    std::vector<std::vector<Rcpp::NumericVector> > quant;
    std::vector<Rcpp::IntegerVector> mr;

    std::vector<scran_markers::SummaryBuffers<double, int> > groupwise;
    initialize_summary_buffers(
        num_groups,
        num_genes,
        groupwise,
        compute_summary_min,
        min,
        compute_summary_mean,
        mean,
        compute_summary_median,
        median,
        compute_summary_max,
        max,
        num_quantiles,
        quant,
        compute_summary_min_rank,
        mr
    );

    scran_markers::summarize_effects(num_genes, num_groups, static_cast<const double*>(effects.begin()), groupwise, opt);

    return format_summary_output(
        num_groups,
        compute_summary_min,
        min,
        compute_summary_mean,
        mean,
        compute_summary_median,
        median,
        compute_summary_max,
        max,
        num_quantiles,
        opt.compute_quantiles,
        quant,
        compute_summary_min_rank,
        mr
    );
}
