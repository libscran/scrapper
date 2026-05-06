#include "config.h"

#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "scran_pca/scran_pca.hpp"

#include "utils_block.h"
#include "utils_other.h"

static Rcpp::NumericMatrix transfer(const Eigen::MatrixXd& x) {
    auto output = create_matrix<Rcpp::NumericMatrix>(x.rows(), x.cols());
    static_assert(!Eigen::MatrixXd::IsRowMajor);
    std::copy_n(x.data(), output.size(), output.begin());
    return output;
}

static Rcpp::NumericVector transfer(const Eigen::VectorXd& x) {
    auto output = sanisizer::create<Rcpp::NumericVector>(x.size());
    std::copy(x.begin(), x.end(), output.begin());
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List run_pca(
    SEXP x,
    Rcpp::RObject number,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject components_from_residuals,
    Rcpp::RObject scale,
    Rcpp::Nullable<Rcpp::IntegerVector> subset,
    Rcpp::RObject realized,
    Rcpp::RObject irlba_work,
    Rcpp::RObject irlba_iterations,
    Rcpp::RObject irlba_tolerance,
    Rcpp::RObject irlba_seed,
    Rcpp::RObject num_threads
) {
    auto mat = Rtatami::BoundNumericPointer(x);
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    irlba::Options iopt;
    set_optional_integer(irlba_work, iopt.extra_work, "extra.work");
    set_integer(irlba_iterations, iopt.max_iterations, "iterations");
    set_number(irlba_tolerance, iopt.convergence_tolerance, "tolerance");
    set_integer(irlba_seed, iopt.seed, "seed");
    iopt.cap_number = true;

    const auto fill_common_options = [&](auto& opt) -> void {
        set_integer(number, opt.number, "number");
        set_bool(scale, opt.scale, "scale");
        set_bool(realized, opt.realize_matrix, "realized");
        opt.irlba_options = iopt;
        set_integer(num_threads, opt.num_threads, "num.threads");
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
            Rcpp::Named("metrics") = Rcpp::List::create(
                Rcpp::Named("converged") = Rcpp::LogicalVector::create(out.metrics.converged),
                Rcpp::Named("iterations") = Rcpp::LogicalVector::create(out.metrics.iterations),
                Rcpp::Named("multiplications") = Rcpp::LogicalVector::create(out.metrics.multiplications)
            )
        );
    };

    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), mat->ptr->ncol())) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        const auto fill_block_options = [&](auto& opt) -> void {
            fill_common_options(opt);
            set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
            set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");
            set_bool(components_from_residuals, opt.components_from_residuals, "components.from.residuals");
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

//[[Rcpp::export(rng=false)]]
Rcpp::List run_pca_defaults() {
    Rcpp::List output;

    {
        scran_pca::SimplePcaOptions opt;
        output["number"] = opt.number;
        output["scale"] = opt.scale;

        if (opt.irlba_options.extra_work.has_value()) {
            output["extra.work"] = R_NilValue;
        } else {
            throw std::runtime_error("unexpected extra.work default for runPca"); 
        }

        output["iterations"] = opt.irlba_options.max_iterations;
        output["tolerance"] = opt.irlba_options.convergence_tolerance;
        output["seed"] = opt.irlba_options.seed;
        output["realized"] = opt.realize_matrix;
        output["num.threads"] = opt.num_threads;
    }

    {
        scran_pca::BlockedPcaOptions opt;
        report_block_weight_policy_default(output, opt.block_weight_policy, "block.weight.policy", "runPca");
        report_variable_block_weight_default(output, opt.variable_block_weight_parameters, "variable.block.weight");
        output["components.from.residuals"] = opt.components_from_residuals;
    }

    return output;
}

