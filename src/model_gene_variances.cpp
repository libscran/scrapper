//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "Rtatami.h"
#include "scran_variances/scran_variances.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List model_gene_variances(
    SEXP x,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    size_t nblocks,
    std::string block_average_policy,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double block_quantile,
    bool mean_filter,
    double min_mean,
    bool transform,
    double span,
    bool use_min_width,
    double min_width,
    int min_window_count,
    int num_threads)
{
    scran_variances::ModelGeneVariancesOptions opt;
    opt.fit_variance_trend_options.mean_filter = mean_filter;
    opt.fit_variance_trend_options.minimum_mean = min_mean;
    opt.fit_variance_trend_options.transform = transform;
    opt.fit_variance_trend_options.span = span;
    opt.fit_variance_trend_options.use_minimum_width = use_min_width;
    opt.fit_variance_trend_options.minimum_width = min_width;
    opt.fit_variance_trend_options.minimum_window_count = min_window_count;
    opt.num_threads = num_threads;

    if (block_average_policy == "mean") {
        opt.block_average_policy = scran_variances::BlockAveragePolicy::MEAN;
    } else if (block_average_policy == "quantile") {
        opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
    } else {
        throw std::runtime_error("block average policy should be either 'mean' or 'quantile'");
    }

    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.block_quantile = block_quantile;

    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    size_t nc = mat->ncol();
    size_t nr = mat->nrow();

    Rcpp::NumericVector means(nr), variances(nr), fitted(nr), residuals(nr);
    scran_variances::ModelGeneVariancesBuffers<double> buffers;
    buffers.means = means.begin();
    buffers.variances = variances.begin();
    buffers.fitted = fitted.begin();
    buffers.residuals = residuals.begin();

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (block_info.size() != nc) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_variances::ModelGeneVariancesBlockedBuffers<double> bbuffers;
        bbuffers.average = buffers;
        bbuffers.per_block.resize(nblocks);

        std::vector<Rcpp::NumericVector> block_mean, block_var, block_fit, block_res;
        block_mean.reserve(nblocks);
        block_var.reserve(nblocks);
        block_fit.reserve(nblocks);
        block_res.reserve(nblocks);

        for (size_t b = 0; b < nblocks; ++b) {
            block_mean.emplace_back(nr);
            bbuffers.per_block[b].means = block_mean.back().begin();
            block_var.emplace_back(nr);
            bbuffers.per_block[b].variances = block_var.back().begin();
            block_fit.emplace_back(nr);
            bbuffers.per_block[b].fitted = block_fit.back().begin();
            block_res.emplace_back(nr);
            bbuffers.per_block[b].residuals = block_res.back().begin();
        }

        scran_variances::model_gene_variances_blocked(*mat, ptr, bbuffers, opt);

        Rcpp::List pb(nblocks);
        for (size_t b = 0; b < nblocks; ++b) {
            pb[b] = Rcpp::List::create(
                Rcpp::Named("means") = block_mean[b],
                Rcpp::Named("variances") = block_var[b],
                Rcpp::Named("fitted") = block_fit[b],
                Rcpp::Named("residuals") = block_res[b]
            );
        }

        return Rcpp::List::create(
            Rcpp::Named("means") = means,
            Rcpp::Named("variances") = variances,
            Rcpp::Named("fitted") = fitted,
            Rcpp::Named("residuals") = residuals,
            Rcpp::Named("per.block") = pb
        );

    } else {
        scran_variances::model_gene_variances(*mat, buffers, opt);
        return Rcpp::List::create(
            Rcpp::Named("means") = means,
            Rcpp::Named("variances") = variances,
            Rcpp::Named("fitted") = fitted,
            Rcpp::Named("residuals") = residuals
        );
    }
}
