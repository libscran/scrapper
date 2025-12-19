#include "config.h"

#include <vector>
#include <stdexcept>
#include <cstddef>

#include "scran_qc/scran_qc.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_qc.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_adt_qc_metrics(SEXP x, Rcpp::List subsets, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto nc = mat->ncol();
    const auto nr = mat->nrow();

    // Setting up the subsets.
    std::vector<Rcpp::LogicalVector> in_subsets;
    std::vector<const int*> sub_ptrs;
    cast_subset_vectors(nr, subsets, in_subsets, sub_ptrs);

    // Creating output containers.
    scran_qc::ComputeAdtQcMetricsBuffers<double, int> buffers;
    Rcpp::NumericVector sums(nc);
    buffers.sum = sums.begin();
    Rcpp::IntegerVector detected(nc);
    buffers.detected = detected.begin();
    std::vector<Rcpp::NumericVector> out_subsets;
    prepare_subset_metrics(nc, sub_ptrs.size(), out_subsets, buffers.subset_sum); 

    // Running QC code.
    scran_qc::ComputeAdtQcMetricsOptions opt;
    opt.num_threads = num_threads;
    scran_qc::compute_adt_qc_metrics(*mat, sub_ptrs, buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("sum") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("subsets") = Rcpp::List(out_subsets.begin(), out_subsets.end())
    );
}

class ConvertedAdtQcMetrics {
public:
    ConvertedAdtQcMetrics(Rcpp::List metrics) {
        check_names(metrics, { "sum", "detected", "subsets" });

        sums = metrics["sum"];
        const auto ncells = sums.size();

        detected = metrics["detected"];
        if (!sanisizer::is_equal(ncells, detected.size())) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }

        Rcpp::List tmp(metrics["subsets"]);
        check_subset_metrics(ncells, tmp, subsets);
    }

private:
    Rcpp::NumericVector sums;
    Rcpp::IntegerVector detected;
    std::vector<Rcpp::NumericVector> subsets;

public:
    auto size() const {
        return sums.size();
    }

    auto num_subsets() const {
        return subsets.size();
    }

    auto to_buffer() const {
        scran_qc::ComputeAdtQcMetricsBuffers<const double, const int> buffers;
        buffers.sum = sums.begin();
        buffers.detected = detected.begin();
        for (auto& s : subsets) {
            buffers.subset_sum.push_back(s.begin());
        }
        return buffers;
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List suggest_adt_qc_thresholds(Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block, double min_detected_drop, double num_mads) {
    ConvertedAdtQcMetrics all_metrics(metrics);
    auto buffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();

    scran_qc::ComputeAdtQcFiltersOptions opt;
    opt.detected_num_mads = num_mads;
    opt.subset_sum_num_mads = num_mads;
    opt.detected_min_drop = min_detected_drop;

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        auto filt = scran_qc::compute_adt_qc_filters_blocked(sanisizer::cast<std::size_t>(ncells), buffers, ptr, opt);
        const auto& dout = filt.get_detected();
        return Rcpp::List::create(
            Rcpp::Named("detected") = Rcpp::NumericVector(dout.begin(), dout.end()),
            Rcpp::Named("subsets") = create_subset_filters(filt.get_subset_sum())
        );
    } else {
        auto filt = scran_qc::compute_adt_qc_filters(sanisizer::cast<std::size_t>(ncells), buffers, opt);
        const auto& ssout = filt.get_subset_sum();
        return Rcpp::List::create(
            Rcpp::Named("detected") = Rcpp::NumericVector::create(filt.get_detected()),
            Rcpp::Named("subsets") = Rcpp::NumericVector(ssout.begin(), ssout.end())
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector filter_adt_qc_metrics(Rcpp::List filters, Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block) {
    ConvertedAdtQcMetrics all_metrics(metrics);
    auto mbuffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();
    const auto nsubs = all_metrics.num_subsets();

    check_names(filters, { "detected", "subsets" });

    Rcpp::LogicalVector keep(ncells);
    auto kptr = static_cast<int*>(keep.begin());

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_qc::AdtQcBlockedFilters filt;

        Rcpp::NumericVector detected(filters["detected"]);
        const auto nblocks = detected.size();
        copy_filters_blocked(nblocks, filters["detected"], filt.get_detected());
        copy_subset_filters_blocked(nsubs, nblocks, filters["subsets"], filt.get_subset_sum());

        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, ptr, kptr);

    } else {
        scran_qc::AdtQcFilters filt;
        filt.get_detected() = parse_filter_unblocked(filters["detected"], "filters$detected");
        copy_subset_filters_unblocked(nsubs, filters["subsets"], filt.get_subset_sum());
        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, kptr);
    }

    return keep;
}
