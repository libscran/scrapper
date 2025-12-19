#include "config.h"

#include <cstddef>

#include "scran_norm/scran_norm.hpp"
#include "sanisizer/sanisizer.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector sanitize_size_factors(Rcpp::NumericVector size_factors, bool handle_zero, bool handle_negative, bool handle_nan, bool handle_infinite) {
    scran_norm::SanitizeSizeFactorsOptions opt;
    opt.handle_zero = (handle_zero ? scran_norm::SanitizeAction::SANITIZE : scran_norm::SanitizeAction::IGNORE);
    opt.handle_negative = (handle_negative ? scran_norm::SanitizeAction::SANITIZE : scran_norm::SanitizeAction::IGNORE);
    opt.handle_infinite = (handle_infinite ? scran_norm::SanitizeAction::SANITIZE : scran_norm::SanitizeAction::IGNORE);
    opt.handle_nan = (handle_nan ? scran_norm::SanitizeAction::SANITIZE : scran_norm::SanitizeAction::IGNORE);
    auto output = Rcpp::clone(size_factors);
    scran_norm::sanitize_size_factors(sanisizer::cast<std::size_t>(output.size()), static_cast<double*>(output.begin()), opt);
    return output; 
}
