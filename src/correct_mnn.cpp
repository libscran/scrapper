#include "config.h"

#include <string>
#include <memory>
#include <stdexcept>
#include <cstddef>

#include "mnncorrect/mnncorrect.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List correct_mnn(
    Rcpp::NumericMatrix x, 
    Rcpp::IntegerVector block, 
    Rcpp::RObject num_neighbors, 
    Rcpp::RObject num_steps, 
    Rcpp::RObject num_threads,
    Rcpp::RObject merge_policy, 
    SEXP builder
) {
    mnncorrect::Options<int, double> opts;
    set_integer(num_neighbors, opts.num_neighbors, "num.neighbors");
    set_integer(num_steps, opts.num_steps, "num.steps");
    set_integer(num_threads, opts.num_threads, "num.threads"); 

    if (!merge_policy.isNULL()) {
        const std::string mp = parse_single_string(merge_policy, "merge.policy");
        if (mp == "input") {
            opts.merge_policy = mnncorrect::MergePolicy::INPUT;
        } else if (mp == "max-variance" || mp == "variance") {
            opts.merge_policy = mnncorrect::MergePolicy::VARIANCE;
        } else if (mp == "max-rss" || mp == "rss") {
            opts.merge_policy = mnncorrect::MergePolicy::RSS;
        } else if (mp == "max-size" || mp == "size") {
            opts.merge_policy = mnncorrect::MergePolicy::SIZE;
        } else {
            throw std::runtime_error("unknown merge policy");
        }
    }

    BiocNeighbors::BuilderPointer ptr(builder);
    opts.builder = std::shared_ptr<BiocNeighbors::Builder>(std::shared_ptr<BiocNeighbors::Builder>{}, &(*ptr)); // aliasing constructor to no-op destruction.

    if (!sanisizer::is_equal(x.ncol(), block.size())) {
        throw std::runtime_error("length of 'block' should equal the number of columns in 'x'");
    }

    auto output = create_matrix<Rcpp::NumericMatrix>(x.nrow(), x.ncol());
    mnncorrect::compute(
        sanisizer::cast<std::size_t>(x.nrow()),
        sanisizer::cast<int>(x.ncol()),
        static_cast<const double*>(x.begin()),
        static_cast<const int*>(block.begin()),
        static_cast<double*>(output.begin()),
        opts
    );

    return Rcpp::List::create(
        Rcpp::Named("corrected") = output
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List correct_mnn_defaults() {
    Rcpp::List output;
    mnncorrect::Options<int, double> opts;
    output["num.neighbors"] = opts.num_neighbors;
    output["num.steps"] = opts.num_steps;
    output["num.threads"] = opts.num_threads;
    if (opts.merge_policy == mnncorrect::MergePolicy::RSS) {
        output["merge.policy"] = "rss";
    } else {
        throw std::runtime_error("unexpected merge.policy default for correctMnn");
    }
    return output;
}
