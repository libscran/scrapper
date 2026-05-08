#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_variances/scran_variances.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List fit_variance_trend(
    Rcpp::NumericVector means,
    Rcpp::NumericVector variances,
    Rcpp::RObject mean_filter,
    Rcpp::RObject min_mean,
    Rcpp::RObject transform,
    Rcpp::RObject span,
    Rcpp::RObject use_min_width,
    Rcpp::RObject min_width,
    Rcpp::RObject min_window_count,
    Rcpp::RObject num_threads
) {
    scran_variances::FitVarianceTrendOptions opt;
    set_bool(mean_filter, opt.mean_filter, "mean.filter");
    set_number(min_mean, opt.minimum_mean, "min.mean");
    set_bool(transform, opt.transform, "transform");
    set_number(span, opt.span, "span");
    set_bool(use_min_width, opt.use_minimum_width, "use.min.width");
    set_number(min_width, opt.minimum_width, "min.width");
    set_integer(min_window_count, opt.minimum_window_count, "min.window.count");
    set_integer(num_threads, opt.num_threads, "num.threads");

    const auto n = means.size();
    if (!sanisizer::is_equal(n, variances.size())) {
        throw std::runtime_error("'means' and 'variances' should have the same length");
    }

    sanisizer::as_size_type<Rcpp::NumericVector>(n);
    Rcpp::NumericVector fitted(n), residuals(n);
    scran_variances::FitVarianceTrendWorkspace<double> wrk; 
    scran_variances::fit_variance_trend(
        sanisizer::cast<std::size_t>(n),
        static_cast<const double*>(means.begin()),
        static_cast<const double*>(variances.begin()),
        static_cast<double*>(fitted.begin()),
        static_cast<double*>(residuals.begin()),
        wrk,
        opt
    );

    return Rcpp::List::create(
        Rcpp::Named("fitted") = fitted,
        Rcpp::Named("residuals") = residuals
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List fit_variance_trend_defaults() {
    Rcpp::List output; 
    scran_variances::FitVarianceTrendOptions opt;
    output["mean.filter"] = opt.mean_filter;
    output["min.mean"] = opt.minimum_mean;
    output["transform"] = opt.transform;
    output["span"] = opt.span;
    output["use.min.width"] = opt.use_minimum_width;
    output["min.width"] = opt.minimum_width;
    output["min.window.count"] = opt.minimum_window_count;
    output["num.threads"] = opt.num_threads;
    return output;
}
