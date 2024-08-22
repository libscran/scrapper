//#include "config.h"

#include "mumosa/mumosa.hpp"
#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector scale_by_neighbors(Rcpp::List distances) {
    size_t nmod = distances.size();
    std::vector<std::pair<double, double> > values;
    values.reserve(nmod);

    for (size_t x = 0; x < nmod; ++x) {
        Rcpp::NumericVector dist(distances[x]);
        values.push_back(mumosa::compute_distance<int, double>(dist.size(), dist.begin()));
    }

    auto output = mumosa::compute_scale<double>(values);
    return Rcpp::NumericVector(output.begin(), output.end());
}
