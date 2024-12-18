#include "Rcpp.h"
#include "kmeans/kmeans.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_kmeans(
    Rcpp::NumericMatrix data,
    int nclusters,
    std::string init_method,
    std::string refine_method,
    bool var_part_optimize_partition,
    double var_part_size_adjustment,
    int lloyd_iterations,
    int hartigan_wong_iterations,
    int hartigan_wong_quick_transfer_iterations,
    bool hartigan_wong_quit_quick_transfer_failure,
    int seed,
    int nthreads)
{
    int ndim = data.nrow();
    size_t nobs = data.ncol();
    auto ptr = static_cast<const double*>(data.begin());

    Rcpp::NumericMatrix centers(ndim, nclusters);
    Rcpp::IntegerVector clusters(nobs);
    auto center_ptr = static_cast<double*>(centers.begin());
    auto cluster_ptr = static_cast<int*>(clusters.begin());

    std::unique_ptr<kmeans::Initialize<> > iptr;
    if (init_method == "random") {
        auto ptr = new kmeans::InitializeRandom;
        ptr->get_options().seed = seed;
        iptr.reset(ptr);
    } else if (init_method == "kmeans++") {
        auto ptr = new kmeans::InitializeKmeanspp;
        ptr->get_options().num_threads = nthreads;
        ptr->get_options().seed = seed;
        iptr.reset(ptr);;
    } else if (init_method == "var-part") {
        auto ptr = new kmeans::InitializeVariancePartition;
        ptr->get_options().optimize_partition = var_part_optimize_partition;
        ptr->get_options().size_adjustment = var_part_size_adjustment;
        iptr.reset(ptr);
    } else {
        throw std::runtime_error("unknown init_method '" + init_method + "'");
    }

    std::unique_ptr<kmeans::Refine<> > rptr;
    if (refine_method == "lloyd") {
        auto ptr = new kmeans::RefineLloyd;
        ptr->get_options().max_iterations = lloyd_iterations;
        ptr->get_options().num_threads = nthreads;
        rptr.reset(ptr);
    } else if (refine_method == "hartigan-wong") {
        auto ptr = new kmeans::RefineHartiganWong;
        ptr->get_options().max_iterations = hartigan_wong_iterations;
        ptr->get_options().max_quick_transfer_iterations = hartigan_wong_quick_transfer_iterations;
        ptr->get_options().quit_on_quick_transfer_convergence_failure = hartigan_wong_quit_quick_transfer_failure;
        ptr->get_options().num_threads = nthreads;
        rptr.reset(ptr);
    }

    auto out = kmeans::compute(kmeans::SimpleMatrix<double, int, int>(ndim, nobs, ptr), *iptr, *rptr, nclusters, center_ptr, cluster_ptr);
    return Rcpp::List::create(
        Rcpp::Named("clusters") = clusters, 
        Rcpp::Named("centers") = centers,
        Rcpp::Named("iterations") = Rcpp::IntegerVector::create(out.iterations),
        Rcpp::Named("status") = Rcpp::IntegerVector::create(out.status) 
    );
}
