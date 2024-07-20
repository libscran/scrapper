//#include "config.h"

#include <vector>

#include "Rcpp.h"
#include "scran_norm/scran_norm.hpp"

#include "ResolvedBatch.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector center_size_factors(Rcpp::NumericVector size_factors, Rcpp::Nullable<Rcpp::IntegerVector> block, bool lowest) {
    auto block_info = ResolvedBatch(block);
    auto ptr = block_info.get();

    scran_norm::CenterSizeFactorsOptions opt;
    opt.block_mode = (lowest ? scran_norm::CenterBlockMode::LOWEST : scran_norm::CenterBlockMode::PER_BLOCK);
    opt.ignore_invalid = true;

    auto output = Rcpp::clone(size_factors);
    if (ptr) {
        scran_norm::center_size_factors_blocked(output.size(), static_cast<double*>(output.begin()), ptr, NULL, opt);
    } else {
        scran_norm::center_size_factors(output.size(), static_cast<double*>(output.begin()), NULL, opt);
    }

    return output; 
}
