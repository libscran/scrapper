//#include "config.h"

#include "Rcpp.h"
#include "mnncorrect/mnncorrect.hpp"
#include "BiocNeighbors.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List correct_mnn(
    Rcpp::NumericMatrix x, 
    Rcpp::IntegerVector block, 
    int num_neighbors, 
    int num_steps, 
    int num_threads,
    std::string merge_policy, 
    SEXP builder)
{
    mnncorrect::Options<int, double> opts;
    opts.num_neighbors = num_neighbors;
    opts.num_steps = num_steps;
    opts.num_threads = num_threads;

    if (merge_policy == "input") {
        opts.merge_policy = mnncorrect::MergePolicy::INPUT;
    } else if (merge_policy == "max-variance" || merge_policy == "variance") {
        opts.merge_policy = mnncorrect::MergePolicy::VARIANCE;
    } else if (merge_policy == "max-rss" || merge_policy == "rss") {
        opts.merge_policy = mnncorrect::MergePolicy::RSS;
    } else if (merge_policy == "max-size" || merge_policy == "size") {
        opts.merge_policy = mnncorrect::MergePolicy::SIZE;
    } else {
        throw std::runtime_error("unknown merge policy");
    }

    BiocNeighbors::BuilderPointer ptr(builder);
    opts.builder = std::shared_ptr<BiocNeighbors::Builder>(std::shared_ptr<BiocNeighbors::Builder>{}, &(*ptr)); // aliasing constructor to no-op destruction.

    if (x.ncol() != block.size()) {
        throw std::runtime_error("length of 'block' should equal the number of columns in 'x'");
    }

    Rcpp::NumericMatrix output(x.nrow(), x.ncol());
    mnncorrect::compute(
        x.nrow(),
        x.ncol(),
        static_cast<const double*>(x.begin()),
        static_cast<const int*>(block.begin()),
        static_cast<double*>(output.begin()),
        opts
    );

    return Rcpp::List::create(
        Rcpp::Named("corrected") = output
    );
}
