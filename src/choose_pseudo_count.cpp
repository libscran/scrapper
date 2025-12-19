#include "config.h"

#include "scran_norm/scran_norm.hpp"

//[[Rcpp::export(rng=false)]]
double choose_pseudo_count(Rcpp::NumericVector size_factors, double quantile, double max_bias, double min_value) {
    scran_norm::ChoosePseudoCountOptions opt;
    opt.quantile = quantile;
    opt.max_bias = max_bias;
    opt.min_value = min_value;
    return scran_norm::choose_pseudo_count(sanisizer::cast<std::size_t>(size_factors.size()), static_cast<const double*>(size_factors.begin()), opt);
}
