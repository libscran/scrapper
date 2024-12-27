//#include "config.h"

#include <vector>
#include <stdexcept>

#include "scran_aggregate/aggregate_across_genes.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "Rcpp.h"
#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_genes(SEXP x, Rcpp::List sets, bool average, int nthreads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    int NC = mat->ncol();

    // Converting the sets into something nice. We need to make explicit copies
    // to ensure that we resolve any ALTREPs that might be present.
    size_t nsets = sets.size();
    std::vector<std::vector<int> > indices;
    indices.reserve(nsets);
    std::vector<std::vector<double> > weights;
    weights.reserve(nsets);
    std::vector<std::tuple<size_t, const int*, const double*> > converted_sets;
    converted_sets.reserve(nsets);

    for (size_t s = 0; s < nsets; ++s) {
        Rcpp::RObject current = sets[s];
        const double* wptr = NULL;

        if (current.sexp_type() == INTSXP) {
            Rcpp::IntegerVector idx(current);
            indices.emplace_back(idx.begin(), idx.end());
            weights.emplace_back(0);

        } else if (current.sexp_type() == VECSXP) {
            Rcpp::List weighted(current);
            if (weighted.size() != 2) {
                throw std::runtime_error("list entries of 'sets' should be of length 2");
            }
            Rcpp::IntegerVector idx(weighted[0]);
            Rcpp::NumericVector wt(weighted[1]);
            if (idx.size() != wt.size()) {
                throw std::runtime_error("list entries of 'sets' should have vectors of equal length");
            }
            indices.emplace_back(idx.begin(), idx.end());
            weights.emplace_back(wt.begin(), wt.end());
            wptr = weights.back().data();

        } else {
            throw std::runtime_error("unsupported type of 'sets' entry");
        }

        for (auto& ii : indices.back()) {
            --ii;
        }
        converted_sets.emplace_back(indices.back().size(), indices.back().data(), wptr);
    }

    // Constructing the outputs.
    scran_aggregate::AggregateAcrossGenesBuffers<double> buffers;
    buffers.sum.reserve(nsets);
    Rcpp::List output(nsets);
    for (size_t s = 0; s < nsets; ++s) {
        Rcpp::NumericVector current(NC);
        output[s] = current;
        buffers.sum.push_back(static_cast<double*>(current.begin()));
    }

    scran_aggregate::AggregateAcrossGenesOptions opt;
    opt.average = average;
    opt.num_threads = nthreads;
    scran_aggregate::aggregate_across_genes(*mat, converted_sets, buffers, opt);

    return output;
}
