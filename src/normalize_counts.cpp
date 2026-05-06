#include "config.h"

#include <vector>
#include <memory>

#include "scran_norm/scran_norm.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP normalize_counts(SEXP x, Rcpp::NumericVector size_factors, Rcpp::RObject log, Rcpp::RObject pseudo_count, Rcpp::RObject log_base, Rcpp::RObject preserve_sparsity) {
    scran_norm::NormalizeCountsOptions opt;
    set_bool(log, opt.log, "log");
    set_number(pseudo_count, opt.pseudo_count, "pseudo.count");
    set_number(log_base, opt.log_base, "log.base");
    set_bool(preserve_sparsity, opt.preserve_sparsity, "preserve.sparsity");

    Rtatami::BoundNumericPointer mat(x);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = scran_norm::normalize_counts(mat->ptr, std::vector<double>(size_factors.begin(), size_factors.end()), opt);
    output->original = mat->original;
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List normalize_counts_defaults() {
    Rcpp::List output;
    scran_norm::NormalizeCountsOptions opt;
    output["log"] = opt.log;
    output["pseudo_count"] = opt.pseudo_count;
    output["log_base"] = opt.log_base;
    output["preserve.sparsity"] = opt.preserve_sparsity;
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_LogNormalizedMatrix(SEXP seed, Rcpp::NumericVector size_factors, double pseudo_count, double log_base) {
    std::vector<double> recip_sf(size_factors.begin(), size_factors.end());
    for (auto& s : recip_sf) {
        s = 1.0/s;
    }

    Rtatami::BoundNumericPointer mat(seed);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        mat->ptr,
        std::make_shared<scran_norm::DelayedLogNormalizeHelper<double, double, int, std::vector<double> > >(std::move(recip_sf), log_base, pseudo_count)
    );
    output->original = mat->original;
    return output;
}
