#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_aggregate/scran_aggregate.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
SEXP aggregate_across_cells(
    SEXP x,
    Rcpp::IntegerVector groups,
    bool compute_sum,
    bool compute_detected,
    bool compute_median,
    int num_threads
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();

    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("length of 'groups' should be equal to the number of columns in 'x'");
    }

    const int* gptr = groups.begin();
    const auto ncombos = tatami_stats::total_groups(gptr, NC);
    scran_aggregate::AggregateAcrossCellsBuffers<double, int, double> buffers;

    Rcpp::NumericMatrix sums;
    if (compute_sum) {
        sums = create_matrix<Rcpp::NumericMatrix>(NR, ncombos);
        buffers.sums.reserve(ncombos);
        double* osum = sums.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.sums.push_back(osum + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::IntegerMatrix detected;
    if (compute_detected) {
        detected = create_matrix<Rcpp::IntegerMatrix>(NR, ncombos);
        buffers.detected.reserve(ncombos);
        int* odet = detected.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.detected.push_back(odet + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::NumericMatrix medians;
    if (compute_median) {
        medians = create_matrix<Rcpp::NumericMatrix>(NR, ncombos);
        buffers.medians.reserve(ncombos);
        double* omedian = medians.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.medians.push_back(omedian + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    scran_aggregate::AggregateAcrossCellsOptions opt;
    opt.num_threads = num_threads;
    scran_aggregate::aggregate_across_cells(*mat, gptr, buffers, opt);

    Rcpp::RObject sums2 = R_NilValue;
    if (compute_sum) {
        sums2 = sums;
    }

    Rcpp::RObject detected2 = R_NilValue;
    if (compute_detected) {
        detected2 = detected;
    }

    Rcpp::RObject medians2 = R_NilValue;
    if (compute_median) {
        medians2 = medians;
    }

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums2,
        Rcpp::Named("detected") = detected2,
        Rcpp::Named("medians") = medians2
    );
}
