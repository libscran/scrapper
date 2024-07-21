//#include "config.h"

#include <vector>

#include "Rcpp.h"
#include "scran_variances/scran_variances.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::List fit_variance_trend(
    Rcpp::NumericVector means,
    Rcpp::NumericVector variances,
    bool mean_filter,
    double min_mean,
    bool transform,
    double span,
    bool use_min_width,
    double min_width,
    int min_window_count,
    int num_threads)
{
    scran_variances::FitVarianceTrendOptions opt;
    opt.mean_filter = mean_filter;
    opt.minimum_mean = min_mean;
    opt.transform = transform;
    opt.span = span;
    opt.use_minimum_width = use_min_width;
    opt.minimum_width = min_width;
    opt.minimum_window_count = min_window_count;
    opt.num_threads = num_threads;

    size_t n = means.size();
    if (n != variances.size()) {
        throw std::runtime_error("'means' and 'variances' should have the same length");
    }

    Rcpp::NumericVector fitted(n), residuals(n);
    scran_variances::FitVarianceTrendWorkspace<double> wrk; 
    scran_variances::fit_variance_trend(
        n,
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
