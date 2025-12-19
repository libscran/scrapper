#include "config.h"

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "scran_qc/scran_qc.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_qc.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_rna_qc_metrics(SEXP x, Rcpp::List subsets, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto nc = mat->ncol();
    const auto nr = mat->nrow();

    // Setting up the subsets.
    std::vector<Rcpp::LogicalVector> in_subsets;
    std::vector<const int*> sub_ptrs;
    cast_subset_vectors(nr, subsets, in_subsets, sub_ptrs);

    // Creating output containers.
    scran_qc::ComputeRnaQcMetricsBuffers<double, int> buffers;

    auto sums = sanisizer::create<Rcpp::NumericVector>(nc);
    buffers.sum = sums.begin();
    auto detected = sanisizer::create<Rcpp::IntegerVector>(nc);
    buffers.detected = detected.begin();
    std::vector<Rcpp::NumericVector> out_subsets;
    prepare_subset_metrics(nc, sub_ptrs.size(), out_subsets, buffers.subset_proportion); 

    // Running QC code.
    scran_qc::ComputeRnaQcMetricsOptions opt;
    opt.num_threads = num_threads;
    scran_qc::compute_rna_qc_metrics(*mat, sub_ptrs, buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("sum") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("subsets") = Rcpp::List(out_subsets.begin(), out_subsets.end())
    );
}

class ConvertedRnaQcMetrics {
public:
    ConvertedRnaQcMetrics(Rcpp::List metrics) {
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
        scran_qc::ComputeRnaQcMetricsBuffers<const double, const int, const double> buffers;
        buffers.sum = sums.begin();
        buffers.detected = detected.begin();
        for (auto& s : subsets) {
            buffers.subset_proportion.push_back(s.begin());
        }
        return buffers;
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List suggest_rna_qc_thresholds(Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block, double num_mads) {
    ConvertedRnaQcMetrics all_metrics(metrics);
    auto buffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();

    scran_qc::ComputeRnaQcFiltersOptions opt;
    opt.sum_num_mads = num_mads;
    opt.detected_num_mads = num_mads;
    opt.subset_proportion_num_mads = num_mads;

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        auto filt = scran_qc::compute_rna_qc_filters_blocked(sanisizer::cast<std::size_t>(ncells), buffers, ptr, opt);
        const auto& sout = filt.get_sum();
        const auto& dout = filt.get_detected();
        return Rcpp::List::create(
            Rcpp::Named("sum") = Rcpp::NumericVector(sout.begin(), sout.end()),
            Rcpp::Named("detected") = Rcpp::NumericVector(dout.begin(), dout.end()),
            Rcpp::Named("subsets") = create_subset_filters(filt.get_subset_proportion())
        );

    } else {
        auto filt = scran_qc::compute_rna_qc_filters(sanisizer::cast<std::size_t>(ncells), buffers, opt);
        const auto& ssout = filt.get_subset_proportion();
        return Rcpp::List::create(
            Rcpp::Named("sum") = Rcpp::NumericVector::create(filt.get_sum()),
            Rcpp::Named("detected") = Rcpp::NumericVector::create(filt.get_detected()),
            Rcpp::Named("subsets") = Rcpp::NumericVector(ssout.begin(), ssout.end())
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector filter_rna_qc_metrics(Rcpp::List filters, Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block) {
    ConvertedRnaQcMetrics all_metrics(metrics);
    auto mbuffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();
    const auto nsubs = all_metrics.num_subsets();

    check_names(metrics, { "sum", "detected", "subsets" });

    auto keep = sanisizer::create<Rcpp::LogicalVector>(ncells);
    auto kptr = static_cast<int*>(keep.begin());

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_qc::RnaQcBlockedFilters filt;

        Rcpp::NumericVector sum(filters["sum"]);
        const auto nblocks = sum.size();
        copy_filters_blocked(nblocks, sum, filt.get_sum());
        copy_filters_blocked(nblocks, filters["detected"], filt.get_detected());
        copy_subset_filters_blocked(nsubs, nblocks, filters["subsets"], filt.get_subset_proportion());

        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, ptr, kptr);

    } else {
        scran_qc::RnaQcFilters filt;
        filt.get_sum() = parse_filter_unblocked(filters["sum"], "filters$sum");
        filt.get_detected() = parse_filter_unblocked(filters["detected"], "filters$detected");
        copy_subset_filters_unblocked(nsubs, filters["subsets"], filt.get_subset_proportion());
        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, kptr);
    }

    return keep;
}
