#include "Rcpp.h"
#include "utils_block.h"
#include "scran_blocks/scran_blocks.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_block_weights(Rcpp::NumericVector sizes, std::string policy, Rcpp::NumericVector variable_block_weight) {
    Rcpp::NumericVector weights(sizes.size());
    scran_blocks::compute_weights(
        sizes.size(),
        static_cast<const double*>(sizes.begin()), 
        parse_block_weight_policy(policy),
        parse_variable_block_weight(variable_block_weight),
        static_cast<double*>(weights.begin())
    );
    return weights;
}
