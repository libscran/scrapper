#include "config.h"

#include <stdexcept>

#include "scran_variances/scran_variances.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector choose_highly_variable_genes(
    Rcpp::NumericVector stats,
    Rcpp::Nullable<Rcpp::IntegerVector> top,
    Rcpp::Nullable<Rcpp::LogicalVector> larger,
    Rcpp::Nullable<Rcpp::LogicalVector> keep_ties,
    Rcpp::Nullable<Rcpp::NumericVector> bound
) {
    scran_variances::ChooseHighlyVariableGenesOptions opt;
    set_integer(top, opt.top, "top");
    set_bool(larger, opt.larger, "larger");
    set_bool(keep_ties, opt.keep_ties, "keep.ties");

    set_number(bound, opt.bound, "bound");
    if (std::isnan(opt.bound)) {
        opt.use_bound = false;
    }

    auto res = scran_variances::choose_highly_variable_genes_index(sanisizer::cast<std::size_t>(stats.size()), static_cast<const double*>(stats.begin()), opt);
    return Rcpp::IntegerVector(res.begin(), res.end());
}

//[[Rcpp::export(rng=false)]]
Rcpp::List choose_highly_variable_genes_defaults() {
    Rcpp::List output;
    scran_variances::ChooseHighlyVariableGenesOptions opt;
    output["top"] = opt.top;
    output["larger"] = opt.larger;
    output["keep.ties"] = opt.keep_ties;
    if (!opt.use_bound) {
        throw std::runtime_error("unexpected bound default for choose_highly_variable_genes");
    }
    output["bound"] = opt.bound;
    return output;
}
