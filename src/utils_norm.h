#ifndef UTILS_NORM_H
#define UTILS_NORM_H

#include "Rcpp.h"

#include "scran_norm/scran_norm.hpp"
#include "utils_other.h"

inline void set_block_mode(const Rcpp::Nullable<Rcpp::CharacterVector>& mode, scran_norm::CenterBlockMode& target) {
    if (mode.isNull()) {
        return;
    }

    const std::string md = parse_single_string(Rcpp::CharacterVector(mode), "mode");
    if (md == "lowest") {
        target = scran_norm::CenterBlockMode::LOWEST;
    } else if (md == "per-block") {
        target = scran_norm::CenterBlockMode::PER_BLOCK;
    } else {
        throw std::runtime_error("unknown value for 'mode'");
    }
}

#endif
