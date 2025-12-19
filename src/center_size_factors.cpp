#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_norm/scran_norm.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector center_size_factors(Rcpp::NumericVector size_factors, Rcpp::Nullable<Rcpp::IntegerVector> block, bool lowest) {
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    scran_norm::CenterSizeFactorsOptions opt;
    opt.block_mode = (lowest ? scran_norm::CenterBlockMode::LOWEST : scran_norm::CenterBlockMode::PER_BLOCK);
    opt.ignore_invalid = true;

    const auto ncells = size_factors.size();
    auto output = Rcpp::clone(size_factors);
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        scran_norm::center_size_factors_blocked(sanisizer::cast<std::size_t>(ncells), static_cast<double*>(output.begin()), ptr, NULL, opt);
    } else {
        scran_norm::center_size_factors(sanisizer::cast<std::size_t>(output.size()), static_cast<double*>(output.begin()), NULL, opt);
    }

    return output; 
}
