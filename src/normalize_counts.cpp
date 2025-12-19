#include "config.h"

#include <vector>

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

