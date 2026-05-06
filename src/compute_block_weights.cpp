#include "config.h"

#include <string>

#include "scran_blocks/scran_blocks.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_block_weights(Rcpp::NumericVector sizes, std::string policy, Rcpp::NumericVector variable_block_weight) {
    auto weights = sanisizer::create<Rcpp::NumericVector>(sizes.size());

    scran_blocks::VariableWeightParameters output;
    if (variable_block_weight.isNULL()) {
        output = parse_variable_block_weight(variable_block_weight, "variable.block.weight");
    }

    scran_blocks::compute_weights(
        sanisizer::cast<std::size_t>(sizes.size()),
        static_cast<const double*>(sizes.begin()), 
        parse_block_weight_policy(policy),
        output,
        static_cast<double*>(weights.begin())
    );

    return weights;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List compute_block_weights_defaults() {
    Rcpp::List output; 
    report_variable_block_weight_default(output, {}, "variable.block.weight");
    return output;
}
