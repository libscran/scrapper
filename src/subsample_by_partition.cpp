#include "config.h"

#include "partisub/partisub.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP subsample_by_partition(Rcpp::IntegerVector partitions, int target, Rcpp::RObject seed, Rcpp::RObject force_non_empty) {
    partisub::Options opt;
    set_integer(seed, opt.seed, "seed");
    set_bool(force_non_empty, opt.force_non_empty, "force.non.empty");

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

//[[Rcpp::export(rng=false)]]
SEXP subsample_by_partition_defaults() {
    Rcpp::List output;
    partisub::Options opt;
    output["seed"] = opt.seed;
    output["force.non.empty"] = opt.force_non_empty;
    return output;
}
