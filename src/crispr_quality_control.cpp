#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_qc/scran_qc.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_qc.h"

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_crispr_qc_metrics(SEXP x, int num_threads) {
    auto raw_mat = Rtatami::BoundNumericPointer(x);
    const auto& mat = raw_mat->ptr;
    const auto nc = mat->ncol();

    // Creating output containers.
    sanisizer::as_size_type<Rcpp::NumericVector>(nc);
    Rcpp::NumericVector sums(nc), max_value(nc);
    sanisizer::as_size_type<Rcpp::IntegerVector>(nc);
    Rcpp::IntegerVector detected(nc), max_index(nc);

    scran_qc::ComputeCrisprQcMetricsBuffers<double, int, double, int> buffers;
    buffers.sum = sums.begin();
    buffers.detected = detected.begin();
    buffers.max_value = max_value.begin();
    buffers.max_index = max_index.begin();

    // Running QC code.
    scran_qc::ComputeCrisprQcMetricsOptions opt;
    opt.num_threads = num_threads;
    scran_qc::compute_crispr_qc_metrics(*mat, buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("sum") = sums,
        Rcpp::Named("detected") = detected,
        Rcpp::Named("max.value") = max_value,
        Rcpp::Named("max.index") = max_index
    );
}

class ConvertedCrisprQcMetrics {
public:
    ConvertedCrisprQcMetrics(Rcpp::List metrics) {
        check_names(metrics, { "sum", "detected", "max.value", "max.index" });

        sums = metrics["sum"];
        const auto ncells = sums.size();

        detected = metrics["detected"];
        if (!sanisizer::is_equal(ncells, detected.size())) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }

        max_value = metrics["max.value"];
        if (!sanisizer::is_equal(ncells, max_value.size())) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }

        max_index = metrics["max.index"];
        if (!sanisizer::is_equal(ncells, max_index.size())) {
            throw std::runtime_error("all 'metrics' vectors should have the same length");
        }
    }

private:
    Rcpp::NumericVector sums, max_value;
    Rcpp::IntegerVector detected, max_index;

public:
    auto size() const {
        return sums.size();
    }

    auto to_buffer() const {
        scran_qc::ComputeCrisprQcMetricsBuffers<const double, const int, const double, const int> buffers;
        buffers.sum = sums.begin();
        buffers.detected = detected.begin();
        buffers.max_value = max_value.begin();
        buffers.max_index = max_index.begin();
        return buffers;
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List suggest_crispr_qc_thresholds(Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block, double num_mads) {
    ConvertedCrisprQcMetrics all_metrics(metrics);
    auto buffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();

    scran_qc::ComputeCrisprQcFiltersOptions opt;
    opt.max_value_num_mads = num_mads;

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        auto filt = scran_qc::compute_crispr_qc_filters_blocked(sanisizer::cast<std::size_t>(ncells), buffers, ptr, opt);
        const auto& mout = filt.get_max_value();
        return Rcpp::List::create(
            Rcpp::Named("max.value") = Rcpp::NumericVector(mout.begin(), mout.end())
        );
    } else {
        auto filt = scran_qc::compute_crispr_qc_filters(sanisizer::cast<std::size_t>(ncells), buffers, opt);
        return Rcpp::List::create(
            Rcpp::Named("max.value") = Rcpp::NumericVector::create(filt.get_max_value())
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector filter_crispr_qc_metrics(Rcpp::List filters, Rcpp::List metrics, Rcpp::Nullable<Rcpp::IntegerVector> block) {
    ConvertedCrisprQcMetrics all_metrics(metrics);
    auto mbuffers = all_metrics.to_buffer();
    const auto ncells = all_metrics.size();

    check_names(metrics, { "max.value" });

    auto keep = sanisizer::create<Rcpp::LogicalVector>(ncells);
    auto kptr = static_cast<int*>(keep.begin());

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_qc::CrisprQcBlockedFilters filt;

        Rcpp::NumericVector max_value(filters["max.value"]);
        const auto nblocks = max_value.size();
        copy_filters_blocked(nblocks, max_value, filt.get_max_value());

        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, ptr, kptr);

    } else {
        scran_qc::CrisprQcFilters filt;
        filt.get_max_value() = parse_filter_unblocked(filters["max.value"], "filters$max.value");
        filt.filter(sanisizer::cast<std::size_t>(ncells), mbuffers, kptr);
    }

    return keep;
}
