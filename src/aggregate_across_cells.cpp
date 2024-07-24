//#include "config.h"

#include <vector>
#include <stdexcept>

#include "scran_aggregate/scran_aggregate.hpp"
#include "tatami_stats/tatami_stats.hpp"

#include "Rcpp.h"
#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_cells(SEXP x, Rcpp::IntegerVector combined, int nthreads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    size_t NC = mat->ncol();
    size_t NR = mat->nrow();

    if (combined.size() != NC) {
        throw std::runtime_error("length of 'combined' should be equal to the number of columns in 'x'");
    }

    size_t ncombos = tatami_stats::total_groups<int>(combined.begin(), NC);
    Rcpp::NumericMatrix sums(NR, ncombos);
    Rcpp::IntegerMatrix detected(NR, ncombos);

    scran_aggregate::AggregateAcrossCellsBuffers<double, int> buffers;
    {
        buffers.sums.reserve(ncombos);
        buffers.detected.reserve(ncombos);
        double* osum = sums.begin();
        int* odet = detected.begin();
        for (size_t i = 0; i < ncombos; ++i, osum += NR, odet += NR) {
            buffers.sums.push_back(osum);
            buffers.detected.push_back(odet);
        }
    }

    scran_aggregate::AggregateAcrossCellsOptions opt;
    opt.num_threads = nthreads;
    scran_aggregate::aggregate_across_cells(*mat, static_cast<const int*>(combined.begin()), buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("detected") = detected
    );
}
