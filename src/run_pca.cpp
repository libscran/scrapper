//#include "config.h"

#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "Rcpp.h"
#include "Rtatami.h"
#include "scran_pca/scran_pca.hpp"

#include "utils_block.h"

static Rcpp::NumericMatrix transfer(const Eigen::MatrixXd& x) {
    Rcpp::NumericMatrix output(x.rows(), x.cols());
    static_assert(!Eigen::MatrixXd::IsRowMajor);
    std::copy_n(x.data(), output.size(), output.begin());
    return output;
}

static Rcpp::NumericVector transfer(const Eigen::VectorXd& x) {
    Rcpp::NumericVector output(x.size());
    std::copy(x.begin(), x.end(), output.begin());
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List run_pca(
    SEXP x,
    int number,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    bool components_from_residuals,
    bool scale,
    Rcpp::Nullable<Rcpp::IntegerVector> subset,
    bool realized,
    int irlba_work,
    int irlba_iterations,
    int irlba_seed,
    int num_threads
) {
    auto mat = Rtatami::BoundNumericPointer(x);
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    irlba::Options iopt;
    iopt.extra_work = irlba_work;
    iopt.max_iterations = irlba_iterations;
    iopt.seed = irlba_seed;
    iopt.cap_number = true;

    const auto fill_common_options = [&](auto& opt) -> void {
        opt.number = number;
        opt.scale = scale;
        opt.realize_matrix = realized;
        opt.irlba_options = iopt;
        opt.num_threads = num_threads;
    };

    Rcpp::List output;
    const auto deposit_outputs = [&](const auto& out) -> Rcpp::List {
        return Rcpp::List::create(
            Rcpp::Named("components") = transfer(out.components),
            Rcpp::Named("rotation") = transfer(out.rotation),
            Rcpp::Named("variance.explained") = transfer(out.variance_explained),
            Rcpp::Named("total.variance") = Rcpp::NumericVector::create(out.total_variance),
            Rcpp::Named("center") = transfer(out.center),
            Rcpp::Named("scale") = transfer(out.scale),
            Rcpp::Named("converged") = Rcpp::LogicalVector::create(out.converged)
        );
    };

    if (ptr) {
        if (block_info.size() != static_cast<size_t>(mat->ptr->ncol())) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        const auto fill_block_options = [&](auto& opt) -> void {
            fill_common_options(opt);
            opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
            opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
            opt.components_from_residuals = components_from_residuals;
        };

        if (subset.isNull()) {
            scran_pca::BlockedPcaOptions opt;
            fill_block_options(opt);
            auto res = scran_pca::blocked_pca(*(mat->ptr), ptr, opt);
            output = deposit_outputs(res);
        } else {
            scran_pca::SubsetPcaBlockedOptions opt;
            fill_block_options(opt);
            auto res = scran_pca::subset_pca_blocked(*(mat->ptr), Rcpp::IntegerVector(subset), ptr, opt);
            output = deposit_outputs(res);
        }

    } else {
        if (subset.isNull()) {
            scran_pca::SimplePcaOptions opt;
            fill_common_options(opt);
            auto res = scran_pca::simple_pca(*(mat->ptr), opt);
            output = deposit_outputs(res);
        } else {
            scran_pca::SubsetPcaOptions opt;
            fill_common_options(opt);
            auto res = scran_pca::subset_pca(*(mat->ptr), Rcpp::IntegerVector(subset), opt);
            output = deposit_outputs(res);
        }
    }

    return output;
}
