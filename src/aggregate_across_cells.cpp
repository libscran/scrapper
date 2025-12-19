#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_aggregate/scran_aggregate.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_cells(SEXP x, Rcpp::IntegerVector groups, int nthreads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();

    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("length of 'groups' should be equal to the number of columns in 'x'");
    }

    const int* gptr = groups.begin();
    const auto ncombos = tatami_stats::total_groups(gptr, NC);

    auto sums = create_matrix<Rcpp::NumericMatrix>(NR, ncombos);
    auto detected = create_matrix<Rcpp::IntegerMatrix>(NR, ncombos);
    scran_aggregate::AggregateAcrossCellsBuffers<double, int> buffers;
    {
        buffers.sums.reserve(ncombos);
        buffers.detected.reserve(ncombos);
        double* osum = sums.begin();
        int* odet = detected.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            const auto offset = sanisizer::product_unsafe<std::size_t>(NR, i);
            buffers.sums.push_back(osum + offset);
            buffers.detected.push_back(odet + offset);
        }
    }

    scran_aggregate::AggregateAcrossCellsOptions opt;
    opt.num_threads = nthreads;
    scran_aggregate::aggregate_across_cells(*mat, gptr, buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums,
        Rcpp::Named("detected") = detected
    );
}
