#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_norm/scran_norm.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_other.h"
#include "utils_norm.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector center_size_factors(Rcpp::NumericVector size_factors, Rcpp::Nullable<Rcpp::IntegerVector> block, int num_blocks, Rcpp::RObject mode) {
    const auto ncells = size_factors.size();
    auto output = Rcpp::clone(size_factors);

    auto block_info = MaybeBlock(block);
    auto block_ptr = block_info.get();
    if (block_ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_norm::CenterSizeFactorsBlockedOptions opt;
        opt.ignore_invalid = true;
        set_block_mode(mode, opt.block_mode);
        scran_norm::center_size_factors_blocked(
            sanisizer::cast<std::size_t>(ncells),
            static_cast<double*>(output.begin()),
            block_ptr,
            static_cast<std::size_t>(num_blocks),
            opt
        );

    } else {
        scran_norm::CenterSizeFactorsOptions opt;
        opt.ignore_invalid = true;
        scran_norm::center_size_factors(sanisizer::cast<std::size_t>(ncells), static_cast<double*>(output.begin()), opt);
    }

    return output; 
}

//[[Rcpp::export(rng=false)]]
Rcpp::List center_size_factors_defaults() {
    Rcpp::List output;
    scran_norm::CenterSizeFactorsBlockedOptions opt;
    if (opt.block_mode == scran_norm::CenterBlockMode::LOWEST) {
        output["mode"] = "lowest";
    } else {
        throw std::runtime_error("unexpected mode default for centerSizeFactors");
    }
    return output;
}
