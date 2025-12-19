#ifndef UNPACK_NEIGHBORS_H
#define UNPACK_NEIGHBORS_H

#include "config.h"

#include <vector>

#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

template<typename Index_ = int, class Distance_ = double>
std::vector<std::vector<std::pair<Index_, Distance_> > > unpack_neighbors(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist) {
    const auto nobs = nnidx.cols();
    const auto nneighbors = nnidx.rows();
    const auto iptr = static_cast<const int*>(nnidx.begin());
    const auto dptr = static_cast<const double*>(nndist.begin());

    auto neighbors = sanisizer::create<std::vector<std::vector<std::pair<Index_, Distance_> > > >(nobs);
    for (I<decltype(nobs)> i = 0; i < nobs; ++i) {
        auto& current = neighbors[i];
        current.reserve(nneighbors);
        for (I<decltype(nneighbors)> k = 0; k < nneighbors; ++k) {
            const auto offset = sanisizer::nd_offset<std::size_t>(k, nneighbors, i);
            current.emplace_back(iptr[offset] - 1, dptr[offset]);
        }
    }

    return neighbors;
}

#endif
