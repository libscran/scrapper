//#include "config.h"

#include <vector>

#include "Rcpp.h"
#include "Rtatami.h"
#include "scran_variances/scran_variances.hpp"

#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List model_gene_variances(
    SEXP x,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    size_t nblocks,
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
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

    if (block_weight_policy == "none") {
        opt.block_weight_policy = scran_blocks::WeightPolicy::NONE;
    } else if (block_weight_policy == "equal") {
        opt.block_weight_policy = scran_blocks::WeightPolicy::EQUAL;
    } else if (block_weight_policy == "variable") {
        opt.block_weight_policy = scran_blocks::WeightPolicy::VARIABLE;
    } else {
        throw std::runtime_error("unknown block weight policy '" + block_weight_policy + "'");
    }

    if (variable_block_weight.size() != 2) {
        throw std::runtime_error("'variable_block_weight' must be a numeric vector of length 2");
    }
    opt.variable_block_weight_parameters.lower_bound = variable_block_weight[0];
    opt.variable_block_weight_parameters.upper_bound = variable_block_weight[1];

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

    auto block_info = ResolvedBatch(block);
    auto ptr = block_info.get();
    if (ptr) {
        scran_variances::ModelGeneVariancesBlockedBuffers<double> bbuffers;
        bbuffers.average = buffers;
        bbuffers.per_block.resize(nblocks);

        std::vector<Rcpp::NumericVector> block_mean, block_var, block_fit, block_res;
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
