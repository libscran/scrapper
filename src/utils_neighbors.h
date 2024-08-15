#ifndef UNPACK_NEIGBHBORS_H
#define UNPACK_NEIGBHBORS_H

#include <vector>
#include "Rcpp.h"

template<typename Index_ = int, class Distance_ = double>
std::vector<std::vector<std::pair<Index_, Distance_> > > unpack_neighbors(Rcpp::IntegerMatrix nnidx, Rcpp::NumericMatrix nndist) {
    size_t nobs = nnidx.cols();
    size_t nneighbors = nnidx.rows();
    const int* iptr = static_cast<const int*>(nnidx.begin());
    const double* dptr = static_cast<const double*>(nndist.begin());

    std::vector<std::vector<std::pair<Index_, Distance_> > > neighbors(nobs);
    for (auto& current : neighbors) {
        current.reserve(nneighbors);
        for (size_t k = 0; k < nneighbors; ++k, ++iptr, ++dptr) {
            current.emplace_back(*iptr - 1, *dptr);
        }
    }

    return neighbors;
}

#endif
