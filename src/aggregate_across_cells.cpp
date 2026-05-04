#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_aggregate/scran_aggregate.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List aggregate_across_cells(
    SEXP x,
    Rcpp::IntegerVector groups,
    Rcpp::Nullable<Rcpp::LogicalVector> compute_sum,
    Rcpp::Nullable<Rcpp::LogicalVector> compute_detected,
    Rcpp::Nullable<Rcpp::LogicalVector> compute_median,
    Rcpp::Nullable<Rcpp::IntegerVector> num_threads
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

    scran_aggregate::AggregateAcrossCellsOptions opt;
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_bool(compute_sum, opt.compute_sums, "compute.sum"); // setting these for convenience + consistency, even though they have no effect if a buffer is supplied. 
    set_bool(compute_detected, opt.compute_detected, "compute.detected");
    set_bool(compute_median, opt.compute_medians, "compute.median");

    Rcpp::NumericMatrix sums;
    if (opt.compute_sums) {
        sums = create_matrix<Rcpp::NumericMatrix>(NR, ncombos);
        sanisizer::reserve(buffers.sums, ncombos);
        double* osum = sums.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.sums.push_back(osum + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::IntegerMatrix detected;
    if (opt.compute_detected) {
        detected = create_matrix<Rcpp::IntegerMatrix>(NR, ncombos);
        sanisizer::reserve(buffers.detected, ncombos);
        int* odet = detected.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.detected.push_back(odet + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::NumericMatrix medians;
    if (opt.compute_medians) {
        medians = create_matrix<Rcpp::NumericMatrix>(NR, ncombos);
        sanisizer::reserve(buffers.medians, ncombos);
        double* omedian = medians.begin();
        for (I<decltype(ncombos)> i = 0; i < ncombos; ++i) {
            buffers.medians.push_back(omedian + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    scran_aggregate::aggregate_across_cells(*mat, gptr, buffers, opt);

    Rcpp::RObject sums2 = R_NilValue;
    if (opt.compute_sums) {
        sums2 = sums;
    }

    Rcpp::RObject detected2 = R_NilValue;
    if (opt.compute_detected) {
        detected2 = detected;
    }

    Rcpp::RObject medians2 = R_NilValue;
    if (opt.compute_medians) {
        medians2 = medians;
    }

    return Rcpp::List::create(
        Rcpp::Named("sums") = sums2,
        Rcpp::Named("detected") = detected2,
        Rcpp::Named("medians") = medians2
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List aggregate_across_cells_defaults() {
    Rcpp::List output;
    scran_aggregate::AggregateAcrossCellsOptions opt;
    output["compute.sum"] = opt.compute_sums;
    output["compute.detected"] = opt.compute_detected;
    output["compute.median"] = opt.compute_medians;
    output["num.threads"] = opt.num_threads;
    return output;
}
