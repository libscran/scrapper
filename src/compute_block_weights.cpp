#include "config.h"

#include <string>

#include "scran_blocks/scran_blocks.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_block_weights(Rcpp::NumericVector sizes, std::string policy, Rcpp::NumericVector variable_block_weight) {
    auto weights = sanisizer::create<Rcpp::NumericVector>(sizes.size());
    scran_blocks::compute_weights(
        sanisizer::cast<std::size_t>(sizes.size()),
        static_cast<const double*>(sizes.begin()), 
        parse_block_weight_policy(policy),
        parse_variable_block_weight(variable_block_weight),
        static_cast<double*>(weights.begin())
    );
    return weights;
}
