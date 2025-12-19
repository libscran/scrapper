#include "config.h"

#include <stdexcept>

#include "phyper/phyper.hpp"
#include "subpar/subpar.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector test_enrichment(Rcpp::IntegerVector overlap, int num_interest, Rcpp::IntegerVector set_sizes, int universe, bool log, int num_threads) {
    const auto nsets = overlap.size();
    if (!sanisizer::is_equal(nsets, set_sizes.size())) {
        throw std::runtime_error("'overlap' and 'set_sizes' should have the same length");
    }

    phyper::Options opt;
    opt.upper_tail = true;
    opt.log = log;

    Rcpp::NumericVector output(nsets);
    double* optr = output.begin(); // avoid Rcpp inside the parallel section.
    const int* olptr = overlap.begin();
    const int* ssptr = set_sizes.begin();

    subpar::parallelize(num_threads, nsets, [&](int, I<decltype(nsets)> start, I<decltype(nsets)> length) {
        for (I<decltype(nsets)> s = start, end = start + length; s < end; ++s) {
            optr[s] = phyper::compute(
                olptr[s],
                ssptr[s],
                universe - ssptr[s],
                num_interest,
                opt
            );
        }
    });

    return output;
}
