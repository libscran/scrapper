#include "config.h"

#include <vector>
#include <memory>

#include "scran_norm/scran_norm.hpp"

//[[Rcpp::export(rng=false)]]
SEXP normalize_counts(SEXP x, Rcpp::NumericVector size_factors, bool log, double pseudo_count, double log_base, bool preserve_sparsity) {
    scran_norm::NormalizeCountsOptions opt;
    opt.log = log;
    opt.pseudo_count = pseudo_count;
    opt.log_base = log_base;
    opt.preserve_sparsity = preserve_sparsity;

    Rtatami::BoundNumericPointer mat(x);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = scran_norm::normalize_counts(mat->ptr, std::vector<double>(size_factors.begin(), size_factors.end()), opt);
    output->original = mat->original;
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_LogNormalizedMatrix(SEXP seed, Rcpp::NumericVector size_factors, double pseudo_count, double log_base) {
    Rtatami::BoundNumericPointer mat(seed);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        mat->ptr,
        std::make_shared<scran_norm::DelayedLogNormalizeHelper<double, double, int, std::vector<double> > >(
            std::vector<double>(size_factors.begin(), size_factors.end()),
            log_base,
            pseudo_count
        )
    );
    output->original = mat->original;
    return output;
}
