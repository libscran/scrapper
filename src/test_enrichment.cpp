#include "Rcpp.h"
#include "phyper/phyper.hpp"
#include "subpar/subpar.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector test_enrichment(Rcpp::IntegerVector overlap, int num_interest, Rcpp::IntegerVector set_sizes, int universe, bool log, int num_threads) {
    size_t nsets = overlap.size();
    if (nsets != static_cast<size_t>(set_sizes.size())) {
        throw std::runtime_error("'overlap' and 'set_sizes' should have the same length");
    }

    phyper::Options opt;
    opt.upper_tail = true;
    opt.log = log;

    Rcpp::NumericVector output(nsets);
    double* optr = output.begin(); // avoid Rcpp inside the parallel section.
    const int* olptr = overlap.begin();
    const int* ssptr = set_sizes.begin();

    subpar::parallelize(num_threads, nsets, [&](int, size_t start, size_t length) {
        for (size_t s = start, end = start + length; s < end; ++s) {
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
