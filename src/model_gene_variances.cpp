#include "config.h"

#include <vector>
#include <stdexcept>
#include <string>

#include "scran_variances/scran_variances.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List model_gene_variances(
    SEXP x,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    int nblocks,
    std::string block_average_policy,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    double block_quantile,
    bool fit_trend,
    bool mean_filter,
    double min_mean,
    bool transform,
    double span,
    bool use_min_width,
    double min_width,
    int min_window_count,
    int num_threads
) {
    scran_variances::ModelGeneVariancesOptions opt;
    opt.fit_variance_trend_options.mean_filter = mean_filter;
    opt.fit_variance_trend_options.minimum_mean = min_mean;
    opt.fit_variance_trend_options.transform = transform;
    opt.fit_variance_trend_options.span = span;
    opt.fit_variance_trend_options.use_minimum_width = use_min_width;
    opt.fit_variance_trend_options.minimum_width = min_width;
    opt.fit_variance_trend_options.minimum_window_count = min_window_count;
    opt.num_threads = num_threads;

    const bool average_block = (block_average_policy != "none");
    if (average_block) {
        if (block_average_policy == "mean") {
            opt.block_average_policy = scran_variances::BlockAveragePolicy::MEAN;
        } else if (block_average_policy == "quantile") {
            opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
        } else {
            throw std::runtime_error("block average policy should be either 'mean' or 'quantile'");
        }
    }

    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.block_quantile = block_quantile;

    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto nc = mat->ncol();
    const auto nr = mat->nrow();

    sanisizer::as_size_type<Rcpp::NumericVector>(nr);
    Rcpp::NumericVector means(average_block ? nr : 0),
        variances(average_block ? nr : 0),
        fitted(average_block && fit_trend ? nr : 0),
        residuals(average_block && fit_trend ? nr : 0);

    scran_variances::ModelGeneVariancesBuffers<double> buffers;
    if (average_block) {
        buffers.means = means.begin();
        buffers.variances = variances.begin();
        if (fit_trend) {
            buffers.fitted = fitted.begin();
            buffers.residuals = residuals.begin();
        } else {
            buffers.fitted = NULL;
            buffers.residuals = NULL;
        }
    } else {
        buffers.means = NULL;
        buffers.variances = NULL;
        buffers.fitted = NULL;
        buffers.residuals = NULL;
    }

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), nc)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_variances::ModelGeneVariancesBlockedBuffers<double> bbuffers;
        bbuffers.average = buffers;
        sanisizer::resize(bbuffers.per_block, nblocks);

        std::vector<Rcpp::NumericVector> block_mean, block_var, block_fit, block_res;
        sanisizer::reserve(block_mean, nblocks);
        sanisizer::reserve(block_var, nblocks);
        if (fit_trend) {
            sanisizer::reserve(block_fit, nblocks);
            sanisizer::reserve(block_res, nblocks);
        }

        for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
            block_mean.emplace_back(nr);
            bbuffers.per_block[b].means = block_mean.back().begin();
            block_var.emplace_back(nr);
            bbuffers.per_block[b].variances = block_var.back().begin();
            if (fit_trend) {
                block_fit.emplace_back(nr);
                bbuffers.per_block[b].fitted = block_fit.back().begin();
                block_res.emplace_back(nr);
                bbuffers.per_block[b].residuals = block_res.back().begin();
            } else {
                bbuffers.per_block[b].fitted = NULL; 
                bbuffers.per_block[b].residuals = NULL;
            }
        }

        scran_variances::model_gene_variances_blocked(*mat, ptr, bbuffers, opt);

        auto pb = sanisizer::create<Rcpp::List>(nblocks);
        for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
            auto block_out = Rcpp::List::create(
                Rcpp::Named("means") = std::move(block_mean[b]),
                Rcpp::Named("variances") = std::move(block_var[b])
            );
            if (fit_trend) {
                block_out["fitted"] = std::move(block_fit[b]);
                block_out["residuals"] = std::move(block_res[b]);
            }
            pb[b] = std::move(block_out);
        }

        Rcpp::List output;
        if (average_block) {
            output["means"] = std::move(means);
            output["variances"] = std::move(variances);
            if (fit_trend) {
                output["fitted"] = fitted;
                output["residuals"] = residuals;
            }
        }
        output["per.block"] = pb;
        return output;

    } else {
        scran_variances::model_gene_variances(*mat, buffers, opt);
        auto output = Rcpp::List::create(Rcpp::Named("means") = means, Rcpp::Named("variances") = variances);
        if (fit_trend) {
            output["fitted"] = fitted;
            output["residuals"] = residuals;
        }
        return output;
    }
}
