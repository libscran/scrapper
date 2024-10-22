//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_qc/scran_qc.hpp"
#include "Rtatami.h"

#include "utils_block.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_adt_qc_metrics(SEXP x, Rcpp::List subsets, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    size_t nc = mat->ncol();
    size_t nr = mat->nrow();

    // Setting up the subsets.
    size_t nsub = subsets.size();
    std::vector<Rcpp::LogicalVector> in_subsets;
    in_subsets.reserve(nsub);
    for (auto sIt = subsets.begin(); sIt != subsets.end(); ++sIt) {
        in_subsets.emplace_back(*sIt);
        if (nr != static_cast<size_t>(in_subsets.back().size())) {
            throw std::runtime_error("each entry of 'subsets' should have the same length as 'nrow(x)'");
        }
    }

    std::vector<const int*> sub_ptrs;
    sub_ptrs.reserve(nsub);
    for (const auto& s : in_subsets) {
        sub_ptrs.push_back(s.begin());
    }

    // Creating output containers.
    Rcpp::NumericVector sums(nc);
    Rcpp::IntegerVector detected(nc);
    std::vector<Rcpp::NumericVector> out_subsets;
    for (size_t s = 0; s < nsub; ++s) {
        out_subsets.emplace_back(nc);
    }

    scran_qc::ComputeAdtQcMetricsBuffers<double, int> buffers;
    buffers.sum = sums.begin();
    buffers.detected = detected.begin();
    for (auto& s : out_subsets) {
        buffers.subset_sum.push_back(s.begin());
    }

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
        if (metrics.size() != 3) {
            throw std::runtime_error("'metrics' should have the same format as the output of 'computeAdtQcMetrics'");
        }

        sums = metrics["sum"];
        size_t ncells = sums.size();

        detected = metrics["detected"];
        if (ncells != static_cast<size_t>(detected.size())) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }

        Rcpp::List tmp(metrics["subsets"]);
        size_t nsubs = tmp.size();
        subsets.reserve(nsubs);
        for (size_t s = 0; s < nsubs; ++s) {
            subsets.emplace_back(tmp[s]);
            if (static_cast<size_t>(subsets.back().size()) != ncells) {
                throw std::runtime_error("all 'metrics' vectors should have the same length");
            }
        }
    }

private:
    Rcpp::NumericVector sums;
    Rcpp::IntegerVector detected;
    std::vector<Rcpp::NumericVector> subsets;

public:
    size_t size() const {
        return sums.size();
    }

    size_t num_subsets() const {
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
    size_t ncells = all_metrics.size();
    size_t nsubs = all_metrics.num_subsets();

    scran_qc::ComputeAdtQcFiltersOptions opt;
    opt.detected_num_mads = num_mads;
    opt.subset_sum_num_mads = num_mads;
    opt.detected_min_drop = min_detected_drop;

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (block_info.size() != ncells) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        auto filt = scran_qc::compute_adt_qc_filters_blocked(ncells, buffers, ptr, opt);
        const auto& ssout = filt.get_subset_sum();
        Rcpp::List subs(nsubs);
        for (size_t s = 0; s < nsubs; ++s) {
            subs[s] = Rcpp::NumericVector(ssout[s].begin(), ssout[s].end());
        }

        const auto& dout = filt.get_detected();
        return Rcpp::List::create(
            Rcpp::Named("detected") = Rcpp::NumericVector(dout.begin(), dout.end()),
            Rcpp::Named("subsets") = subs
        );
    } else {
        auto filt = scran_qc::compute_adt_qc_filters(ncells, buffers, opt);
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
    size_t ncells = all_metrics.size();
    size_t nsubs = all_metrics.num_subsets();

    if (filters.size() != 2) {
        throw std::runtime_error("'filters' should have the same format as the output of 'suggestAdtQcFilters'");
    }

    Rcpp::LogicalVector keep(ncells);

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (block_info.size() != ncells) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_qc::AdtQcBlockedFilters filt;

        Rcpp::NumericVector detected(filters["detected"]);
        size_t nblocks = detected.size();
        auto& df = filt.get_detected();
        df.insert(df.end(), detected.begin(), detected.end());

        Rcpp::List subsets(filters["subsets"]);
        if (static_cast<size_t>(subsets.size()) != nsubs) {
            throw std::runtime_error("'filters$subsets' should have the same length as the number of subsets in 'metrics'");
        }
        auto& ssf = filt.get_subset_sum();
        ssf.reserve(nsubs);
        for (size_t s = 0; s < nsubs; ++s) {
            Rcpp::NumericVector cursub(subsets[s]);
            if (static_cast<size_t>(cursub.size()) != nblocks) {
                throw std::runtime_error("each entry of 'filters$subsets' should have the same length as 'filters$detected'");
            }
            ssf.emplace_back(cursub.begin(), cursub.end());
        }

        if (block_info.number() > nblocks) {
            throw std::runtime_error("'block' contains out-of-range indices");
        }
        filt.filter(ncells, mbuffers, ptr, static_cast<int*>(keep.begin()));

    } else {
        scran_qc::AdtQcFilters filt;

        Rcpp::NumericVector detected(filters["detected"]);
        if (detected.size() != 1) {
            throw std::runtime_error("'filters$detected' should contain a single threshold");
        }
        filt.get_detected() = detected[0];

        Rcpp::NumericVector subsets(filters["subsets"]);
        if (static_cast<size_t>(subsets.size()) != nsubs) {
            throw std::runtime_error("'filters$subsets' should have the same length as the number of subsets in 'metrics'");
        }
        auto& ssf = filt.get_subset_sum();
        ssf.insert(ssf.end(), subsets.begin(), subsets.end());

        filt.filter(ncells, mbuffers, static_cast<int*>(keep.begin()));
    }

    return keep;
}
