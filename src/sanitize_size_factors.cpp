#include "config.h"

#include <cstddef>
#include <stdexcept>
#include <string>

#include "scran_norm/scran_norm.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

void set_sanitize_option(const Rcpp::RObject& val, scran_norm::SanitizeAction& target, const std::string& arg) {
    if (val.isNULL()) {
        return;
    }
    const auto act = parse_single_string(val, arg);
    if (act == "sanitize") {
        target = scran_norm::SanitizeAction::SANITIZE;
    } else if (act == "error") {
        target = scran_norm::SanitizeAction::ERROR;
    } else if (act == "ignore") {
        target = scran_norm::SanitizeAction::IGNORE;
    } else {
        throw std::runtime_error("unknown option '" + act + "' for '" + arg + "'");
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector sanitize_size_factors(
    Rcpp::NumericVector size_factors,
    Rcpp::RObject handle_zero,
    Rcpp::RObject handle_negative,
    Rcpp::RObject handle_nan,
    Rcpp::RObject handle_infinite
) {
    scran_norm::SanitizeSizeFactorsOptions opt;
    set_sanitize_option(handle_zero, opt.handle_zero, "handle.zero");
    set_sanitize_option(handle_negative, opt.handle_negative, "handle.negative");
    set_sanitize_option(handle_infinite, opt.handle_infinite, "handle.negative");
    set_sanitize_option(handle_nan, opt.handle_nan, "handle.negative");
    auto output = Rcpp::clone(size_factors);
    scran_norm::sanitize_size_factors(sanisizer::cast<std::size_t>(output.size()), static_cast<double*>(output.begin()), opt);
    return output; 
}

//[[Rcpp::export(rng=false)]]
Rcpp::List sanitize_size_factors_defaults() {
    Rcpp::List output;
    scran_norm::SanitizeSizeFactorsOptions opt;

    if (opt.handle_zero != scran_norm::SanitizeAction::ERROR) {
        throw std::runtime_error("expected handle.zero default for sanitizeSizeFactors");
    }
    output["handle.zero"] = "error";

    if (opt.handle_negative != scran_norm::SanitizeAction::ERROR) {
        throw std::runtime_error("expected handle.negative default for sanitizeSizeFactors");
    }
    output["handle.negative"] = "error";

    if (opt.handle_infinite != scran_norm::SanitizeAction::ERROR) {
        throw std::runtime_error("expected handle.infinite default for sanitizeSizeFactors");
    }
    output["handle.infinite"] = "error";

    if (opt.handle_nan != scran_norm::SanitizeAction::ERROR) {
        throw std::runtime_error("expected handle.nan default for sanitizeSizeFactors");
    }
    output["handle.nan"] = "error";

    return output;
}
