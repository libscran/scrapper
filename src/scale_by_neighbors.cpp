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
        std::vector<double> copy(dist.begin(), dist.end()); // creating a copy as compute_distance() will shuffle the distances to compute the median.
        values.push_back(mumosa::compute_distance<int, double>(copy.size(), copy.data()));
    }

    auto output = mumosa::compute_scale<double>(values);
    return Rcpp::NumericVector(output.begin(), output.end());
}
