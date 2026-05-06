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
    Rcpp::RObject block_average_policy,
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject block_quantile,
    Rcpp::RObject fit_trend,
    Rcpp::RObject mean_filter,
    Rcpp::RObject min_mean,
    Rcpp::RObject transform,
    Rcpp::RObject span,
    Rcpp::RObject use_min_width,
    Rcpp::RObject min_width,
    Rcpp::RObject min_window_count,
    Rcpp::RObject num_threads
) {
    scran_variances::ModelGeneVariancesOptions opt;
    set_bool(mean_filter, opt.fit_variance_trend_options.mean_filter, "mean.filter");
    set_number(min_mean, opt.fit_variance_trend_options.minimum_mean, "min.mean");
    set_bool(transform, opt.fit_variance_trend_options.transform, "transform");
    set_number(span, opt.fit_variance_trend_options.span, "span");
    set_bool(use_min_width, opt.fit_variance_trend_options.use_minimum_width, "use.min.width");
    set_number(min_width, opt.fit_variance_trend_options.minimum_width, "min.width");
    set_integer(min_window_count, opt.fit_variance_trend_options.minimum_window_count, "min.window.count");
    set_integer(num_threads, opt.num_threads, "num.threads");

    if (!block_average_policy.isNULL()) {
        const std::string bap = parse_single_string(block_average_policy, "block.average.policy");
        if (bap == "mean") {
            opt.block_average_policy = scran_variances::BlockAveragePolicy::MEAN;
        } else if (bap == "quantile") {
            opt.block_average_policy = scran_variances::BlockAveragePolicy::QUANTILE;
        } else if (bap == "none") {
            opt.block_average_policy = scran_variances::BlockAveragePolicy::NONE;
        } else {
            throw std::runtime_error("block average policy should be either 'mean', 'quantile' or 'none'");
        }
    }

    set_number(block_quantile, opt.block_quantile, "block.quantile");
    set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
    set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");

    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto nc = mat->ncol();
    const auto nr = mat->nrow();
    sanisizer::as_size_type<Rcpp::NumericVector>(nr); // make sure all downstream allocations are safe.

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), nc)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_variances::ModelGeneVariancesBlockedBuffers<double> buffers;

        const bool average_block = (opt.block_average_policy != scran_variances::BlockAveragePolicy::NONE);
        Rcpp::NumericVector means(average_block ? nr : 0),
            variances(average_block ? nr : 0),
            fitted(average_block && fit_trend ? nr : 0),
            residuals(average_block && fit_trend ? nr : 0);

        if (average_block) {
            buffers.average.means = means.begin();
            buffers.average.variances = variances.begin();
            if (fit_trend) {
                buffers.average.fitted = fitted.begin();
                buffers.average.residuals = residuals.begin();
            } else {
                buffers.average.fitted = NULL;
                buffers.average.residuals = NULL;
            }
        } else {
            buffers.average.means = NULL;
            buffers.average.variances = NULL;
            buffers.average.fitted = NULL;
            buffers.average.residuals = NULL;
        }

        sanisizer::resize(buffers.per_block, nblocks);
        std::vector<Rcpp::NumericVector> block_mean, block_var, block_fit, block_res;
        sanisizer::reserve(block_mean, nblocks);
        sanisizer::reserve(block_var, nblocks);
        if (fit_trend) {
            sanisizer::reserve(block_fit, nblocks);
            sanisizer::reserve(block_res, nblocks);
        }

        for (I<decltype(nblocks)> b = 0; b < nblocks; ++b) {
            block_mean.emplace_back(nr);
            buffers.per_block[b].means = block_mean.back().begin();
            block_var.emplace_back(nr);
            buffers.per_block[b].variances = block_var.back().begin();
            if (fit_trend) {
                block_fit.emplace_back(nr);
                buffers.per_block[b].fitted = block_fit.back().begin();
                block_res.emplace_back(nr);
                buffers.per_block[b].residuals = block_res.back().begin();
            } else {
                buffers.per_block[b].fitted = NULL; 
                buffers.per_block[b].residuals = NULL;
            }
        }

        scran_variances::model_gene_variances_blocked(*mat, ptr, buffers, opt);

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
        scran_variances::ModelGeneVariancesBuffers<double> buffers;

        Rcpp::NumericVector means(nr), variances(nr), fitted(fit_trend ? nr : 0), residuals(fit_trend ? nr : 0);
        buffers.means = means.begin();
        buffers.variances = variances.begin();
        if (fit_trend) {
            buffers.fitted = fitted.begin();
            buffers.residuals = residuals.begin();
        } else {
            buffers.fitted = NULL;
            buffers.residuals = NULL;
        }

        scran_variances::model_gene_variances(*mat, buffers, opt);
        auto output = Rcpp::List::create(Rcpp::Named("means") = means, Rcpp::Named("variances") = variances);
        if (fit_trend) {
            output["fitted"] = fitted;
            output["residuals"] = residuals;
        }
        return output;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List model_gene_variances_defaults() {
    Rcpp::List output;
    scran_variances::ModelGeneVariancesOptions opt;
    output["mean.filter"] = opt.fit_variance_trend_options.mean_filter;
    output["min.mean"] = opt.fit_variance_trend_options.minimum_mean;
    output["transform"] = opt.fit_variance_trend_options.transform;
    output["span"] = opt.fit_variance_trend_options.span;
    output["use.min.width"] = opt.fit_variance_trend_options.use_minimum_width;
    output["min.width"] = opt.fit_variance_trend_options.minimum_width;
    output["min.window.count"] = opt.fit_variance_trend_options.minimum_window_count;
    output["num.threads"] = opt.num_threads;

    if (opt.block_average_policy == scran_variances::BlockAveragePolicy::MEAN) {
        output["block.average.policy"] = "mean";
    } else {
        throw std::runtime_error("unexpected block.average.policy default for modelGeneVariances");
    }
    output["block.quantile"] = opt.block_quantile;

    report_block_weight_policy_default(output, opt.block_weight_policy, "block.weight.policy", "modelGeneVariances");
    report_variable_block_weight_default(output, opt.variable_block_weight_parameters, "variable.block.weight");
    return output;
}
