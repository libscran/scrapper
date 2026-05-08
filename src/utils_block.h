#ifndef UTILS_BLOCK_H
#define UTILS_BLOCK_H

#include "config.h"

#include <stdexcept>
#include <string>

#include "scran_blocks/scran_blocks.hpp"

#include "utils_other.h"

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

inline void set_block_weight_policy(const Rcpp::RObject& block_weight_policy, scran_blocks::WeightPolicy& target, const std::string& arg) {
    if (block_weight_policy.isNULL()) {
        return;
    }
    const std::string bwp = parse_single_string(block_weight_policy, arg);
    target = parse_block_weight_policy(bwp);
}

inline scran_blocks::VariableWeightParameters parse_variable_block_weight(const Rcpp::RObject& variable_block_weight, const std::string& arg) {
    if (variable_block_weight.sexp_type() != REALSXP) {
        throw std::runtime_error("'" + arg + "' must be a numeric vector");
    }

    Rcpp::NumericVector vbw(variable_block_weight);
    if (vbw.size() != 2) {
        throw std::runtime_error("'" + arg + "' must be a numeric vector of length 2");
    }

    scran_blocks::VariableWeightParameters output;
    output.lower_bound = vbw[0];
    output.upper_bound = vbw[1];
    return output;
}

inline void set_variable_block_weight(const Rcpp::RObject& variable_block_weight, scran_blocks::VariableWeightParameters& target, const std::string& arg) {
    if (variable_block_weight.isNULL()) {
        return;
    }
    target = parse_variable_block_weight(variable_block_weight, arg);
}

inline void report_block_weight_policy_default(Rcpp::List& output, const scran_blocks::WeightPolicy& def, const std::string& arg, const std::string& fun) {
    if (def == scran_blocks::WeightPolicy::VARIABLE) {
        output["block.weight.policy"] = "variable";
    } else {
        throw std::runtime_error("unexpected " + arg + " default for " + fun);
    }
}

inline void report_variable_block_weight_default(Rcpp::List& output, const scran_blocks::VariableWeightParameters& def, const std::string& arg) {
    output[arg] = Rcpp::NumericVector::create(def.lower_bound, def.upper_bound);
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
