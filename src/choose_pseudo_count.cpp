#include "config.h"

#include "scran_norm/scran_norm.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
double choose_pseudo_count(
    Rcpp::NumericVector size_factors,
    Rcpp::Nullable<Rcpp::NumericVector> quantile,
    Rcpp::Nullable<Rcpp::NumericVector> max_bias,
    Rcpp::Nullable<Rcpp::NumericVector> min_value
) {
    scran_norm::ChoosePseudoCountOptions opt;
    set_number(quantile, opt.quantile, "quantile");
    set_number(max_bias, opt.max_bias, "max.bias");
    set_number(min_value, opt.min_value, "min.value");
    return scran_norm::choose_pseudo_count(sanisizer::cast<std::size_t>(size_factors.size()), static_cast<const double*>(size_factors.begin()), opt);
}

//[[Rcpp::export(rng=false)]]
Rcpp::List choose_pseudo_count_defaults() {
    Rcpp::List output;
    scran_norm::ChoosePseudoCountOptions opt;
    output["quantile"] = opt.quantile;
    output["max.bias"] = opt.max_bias;
    output["min.value"] = opt.min_value;
    return output;
}
