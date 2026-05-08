#include "config.h"

#include "clrm1/clrm1.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_clrm1_factors(SEXP x, Rcpp::RObject num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;

    clrm1::Options opt;
    set_integer(num_threads, opt.num_threads, "num.threads");
    auto output = sanisizer::create<Rcpp::NumericVector>(mat->ncol());
    clrm1::compute(*mat, opt, static_cast<double*>(output.begin()));
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List compute_clrm1_factors_defaults() {
    Rcpp::List output;
    clrm1::Options opt;
    output["num.threads"] = opt.num_threads;
    return output;
}
