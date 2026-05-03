#include "config.h"

#include "partisub/partisub.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP subsample_by_partition(Rcpp::IntegerVector partitions, int target, double seed, bool force_non_empty) {
    partisub::Options opt;
    opt.seed = sanisizer::from_float<unsigned long long>(seed);
    opt.force_non_empty = force_non_empty;

    auto res = partisub::compute(
        sanisizer::cast<int>(partitions.size()), 
        static_cast<const int*>(partitions.begin()),
        target,
        opt
    );

    for (auto& r : res) {
        ++r;
    }

    return Rcpp::IntegerVector(res.begin(), res.end());
}
