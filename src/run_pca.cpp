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
    bool realized,
    int irlba_work,
    int irlba_iterations,
    int irlba_seed,
    int num_threads)
{
    auto mat = Rtatami::BoundNumericPointer(x);
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    irlba::Options iopt;
    iopt.extra_work = irlba_work;
    iopt.max_iterations = irlba_iterations;
    iopt.seed = irlba_seed;

    Rcpp::List output;

    if (ptr) {
        if (block_info.size() != static_cast<size_t>(mat->ptr->ncol())) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_pca::BlockedPcaOptions opt;
        opt.number = number;
        opt.scale = scale;
        opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
        opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
        opt.components_from_residuals = components_from_residuals;
        opt.realize_matrix = realized;
        opt.irlba_options = iopt;
        opt.num_threads = num_threads;

        auto out = scran_pca::blocked_pca(*(mat->ptr), ptr, opt);
        output = Rcpp::List::create(
            Rcpp::Named("components") = transfer(out.components),
            Rcpp::Named("rotation") = transfer(out.rotation),
            Rcpp::Named("variance.explained") = transfer(out.variance_explained),
            Rcpp::Named("total.variance") = Rcpp::NumericVector::create(out.total_variance),
            Rcpp::Named("center") = transfer(out.center),
            Rcpp::Named("scale") = transfer(out.scale)//,
//            Rcpp::Named("converged") = Rcpp::LogicalVector::create(out.scale)
        );

    } else {
        scran_pca::SimplePcaOptions opt;
        opt.number = number;
        opt.scale = scale;
        opt.realize_matrix = realized;
        opt.irlba_options = iopt;
        opt.num_threads = num_threads;

        auto out = scran_pca::simple_pca(*(mat->ptr), opt);
        output = Rcpp::List::create(
            Rcpp::Named("components") = transfer(out.components),
            Rcpp::Named("rotation") = transfer(out.rotation),
            Rcpp::Named("variance.explained") = transfer(out.variance_explained),
            Rcpp::Named("total.variance") = Rcpp::NumericVector::create(out.total_variance),
            Rcpp::Named("center") = transfer(out.center),
            Rcpp::Named("scale") = transfer(out.scale)//,
//            Rcpp::Named("converged") = Rcpp::LogicalVector::create(out.scale)
        );
    }

    return output;
}
