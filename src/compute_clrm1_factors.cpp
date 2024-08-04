//#include "config.h"

#include "clrm1/clrm1.hpp"
#include "Rcpp.h"
#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_clrm1_factors(SEXP x, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;

    clrm1::Options opt;
    opt.num_threads = num_threads;
    Rcpp::NumericVector output(mat->ncol());
    clrm1::compute(*mat, opt, static_cast<double*>(output.begin()));
    return output;
}
