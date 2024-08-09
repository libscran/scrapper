//#include "config.h"

#include "Rcpp.h"
#include "mnncorrect/mnncorrect.hpp"
#include "BiocNeighbors.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List correct_mnn(
    Rcpp::NumericMatrix x, 
    Rcpp::IntegerVector block, 
    int num_neighbors, 
    double num_mads, 
    int robust_iterations,
    double robust_trim,
    int num_threads,
    int mass_cap,
    Rcpp::Nullable<Rcpp::IntegerVector> order, 
    std::string ref_policy, 
    SEXP builder)
{
    mnncorrect::Options<int, int, double> opts;
    opts.num_neighbors = num_neighbors;
    opts.num_mads = num_mads;
    opts.robust_iterations = robust_iterations;
    opts.robust_trim = robust_trim;
    opts.mass_cap = mass_cap;
    opts.num_threads = num_threads;

    if (order.isNotNull()) {
        Rcpp::IntegerVector order_(order);
        opts.order.insert(opts.order.end(), order_.begin(), order_.end());
    }

    if (ref_policy == "input") {
        opts.reference_policy = mnncorrect::ReferencePolicy::INPUT;
    } else if (ref_policy == "max-variance") {
        opts.reference_policy = mnncorrect::ReferencePolicy::MAX_VARIANCE;
    } else if (ref_policy == "max-rss") {
        opts.reference_policy = mnncorrect::ReferencePolicy::MAX_RSS;
    } else if (ref_policy == "max-size") {
        opts.reference_policy = mnncorrect::ReferencePolicy::MAX_SIZE;
    } else {
        throw std::runtime_error("unknown reference policy");
    }

    BiocNeighbors::BuilderPointer ptr(builder);
    opts.builder = std::shared_ptr<BiocNeighbors::Builder>(std::shared_ptr<BiocNeighbors::Builder>{}, &(*ptr)); // aliasing constructor to no-op destruction.

    if (x.ncol() != block.size()) {
        throw std::runtime_error("length of 'block' should equal the number of columns in 'x'");
    }

    Rcpp::NumericMatrix output(x.nrow(), x.ncol());
    auto res = mnncorrect::compute(x.nrow(), x.ncol(), static_cast<const double*>(x.begin()), static_cast<const int*>(block.begin()), static_cast<double*>(output.begin()), opts);

    return Rcpp::List::create(
        Rcpp::Named("corrected") = output,
        Rcpp::Named("merge.order") = Rcpp::IntegerVector(res.merge_order.begin(), res.merge_order.end()),
        Rcpp::Named("num.pairs") = Rcpp::IntegerVector(res.num_pairs.begin(), res.num_pairs.end())
    );
}
