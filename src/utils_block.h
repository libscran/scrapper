#ifndef UTILS_BLOCK_H
#define UTILS_BLOCK_H

#include "config.h"

#include <stdexcept>
#include <string>

#include "scran_blocks/scran_blocks.hpp"

inline scran_blocks::WeightPolicy parse_block_weight_policy(const std::string& block_weight_policy) {
    scran_blocks::WeightPolicy output  = scran_blocks::WeightPolicy::SIZE;
    if (block_weight_policy == "none" || block_weight_policy == "size") {
        ;
    } else if (block_weight_policy == "equal") {
        output = scran_blocks::WeightPolicy::EQUAL;
    } else if (block_weight_policy == "variable") {
        output = scran_blocks::WeightPolicy::VARIABLE;
    } else {
        throw std::runtime_error("unknown block weight policy '" + block_weight_policy + "'");
    }
    return output;
}

inline scran_blocks::VariableWeightParameters parse_variable_block_weight(const Rcpp::NumericVector& variable_block_weight) {
    if (variable_block_weight.size() != 2) {
        throw std::runtime_error("'variable_block_weight' must be a numeric vector of length 2");
    }
    scran_blocks::VariableWeightParameters output;
    output.lower_bound = variable_block_weight[0];
    output.upper_bound = variable_block_weight[1];
    return output;
}

class MaybeBlock {
public:
    MaybeBlock(Rcpp::Nullable<Rcpp::IntegerVector> block) {
        my_has_block = block.isNotNull();
        if (my_has_block) {
            my_block = Rcpp::IntegerVector(block);
        }
        return;
    }

    const int* get() const {
        return (my_has_block ? static_cast<const int*>(my_block.begin()) : NULL);
    }

    auto size() const {
        return my_block.size();
    }

private:
    bool my_has_block;
    // Need to carry along the vector to avoid garbage 
    // collection of any allocated memory (e.g., ALTREP'd)
    // from the Rcpp::Integer initialization.
    Rcpp::IntegerVector my_block;
};


#endif
