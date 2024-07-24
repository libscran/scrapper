//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_norm/scran_norm.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector center_size_factors(Rcpp::NumericVector size_factors, Rcpp::Nullable<Rcpp::IntegerVector> block, bool lowest) {
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    scran_norm::CenterSizeFactorsOptions opt;
    opt.block_mode = (lowest ? scran_norm::CenterBlockMode::LOWEST : scran_norm::CenterBlockMode::PER_BLOCK);
    opt.ignore_invalid = true;

    size_t ncells = size_factors.size();
    auto output = Rcpp::clone(size_factors);
    if (ptr) {
        if (block_info.size() != ncells) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        scran_norm::center_size_factors_blocked(ncells, static_cast<double*>(output.begin()), ptr, NULL, opt);
    } else {
        scran_norm::center_size_factors(output.size(), static_cast<double*>(output.begin()), NULL, opt);
    }

    return output; 
}
