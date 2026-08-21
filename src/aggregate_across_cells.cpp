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
    int num_groups,
    Rcpp::RObject compute_sum,
    Rcpp::RObject compute_detected,
    Rcpp::RObject compute_median,
    Rcpp::RObject num_threads
) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto NC = mat->ncol();
    const auto NR = mat->nrow();

    if (!sanisizer::is_equal(groups.size(), NC)) {
        throw std::runtime_error("length of 'groups' should be equal to the number of columns in 'x'");
    }

    scran_aggregate::AggregateAcrossCellsBuffers<double, int, double> buffers;

    scran_aggregate::AggregateAcrossCellsOptions opt;
    set_integer(num_threads, opt.num_threads, "num.threads");
    set_bool(compute_sum, opt.compute_sum, "compute.sum"); // setting these for convenience + consistency, even though they have no effect if a buffer is supplied. 
    set_bool(compute_detected, opt.compute_detected, "compute.detected");
    set_bool(compute_median, opt.compute_median, "compute.median");

    Rcpp::NumericMatrix sums;
    if (opt.compute_sum) {
        sums = create_matrix<Rcpp::NumericMatrix>(NR, num_groups);
        buffers.sum.reserve(num_groups);
        double* osum = sums.begin();
        for (int i = 0; i < num_groups; ++i) {
            buffers.sum.push_back(osum + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::IntegerMatrix detected;
    if (opt.compute_detected) {
        detected = create_matrix<Rcpp::IntegerMatrix>(NR, num_groups);
        buffers.detected.reserve(num_groups);
        int* odet = detected.begin();
        for (int i = 0; i < num_groups; ++i) {
            buffers.detected.push_back(odet + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    Rcpp::NumericMatrix medians;
    if (opt.compute_median) {
        medians = create_matrix<Rcpp::NumericMatrix>(NR, num_groups);
        buffers.median.reserve(num_groups);
        double* omedian = medians.begin();
        for (int i = 0; i < num_groups; ++i) {
            buffers.median.push_back(omedian + sanisizer::product_unsafe<std::size_t>(NR, i));
        }
    }

    scran_aggregate::aggregate_across_cells(
        *mat,
        static_cast<const int*>(groups.begin()),
        sanisizer::cast<std::size_t>(num_groups),
        buffers,
        opt
    );

    Rcpp::RObject sums2 = R_NilValue;
    if (opt.compute_sum) {
        sums2 = sums;
    }

    Rcpp::RObject detected2 = R_NilValue;
    if (opt.compute_detected) {
        detected2 = detected;
    }

    Rcpp::RObject medians2 = R_NilValue;
    if (opt.compute_median) {
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
    output["compute.sum"] = opt.compute_sum;
    output["compute.detected"] = opt.compute_detected;
    output["compute.median"] = opt.compute_median;
    output["num.threads"] = opt.num_threads;
    return output;
}
