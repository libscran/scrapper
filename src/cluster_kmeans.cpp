#include "config.h"

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

#include "kmeans/kmeans.hpp"
#include "kmeans_tatami/kmeans_tatami.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_kmeans(
    SEXP x,
    int nclusters,
    bool tatami,
    std::string init_method,
    std::string refine_method,
    Rcpp::RObject random_seed,
    Rcpp::RObject kmeanspp_seed,
    Rcpp::RObject var_part_optimize_partition,
    Rcpp::RObject var_part_size_adjustment,
    Rcpp::RObject lloyd_iterations,
    Rcpp::RObject hartigan_wong_iterations,
    Rcpp::RObject hartigan_wong_quick_transfer_iterations,
    Rcpp::RObject hartigan_wong_quit_quick_transfer_failure,
    Rcpp::RObject nthreads
) {
    std::optional<Rcpp::NumericMatrix> data;
    std::unique_ptr<kmeans::Matrix<int, double> > mptr;

    if (!tatami) {
        data = x;
        mptr.reset(new kmeans::SimpleMatrix<int, double>(
            sanisizer::cast<std::size_t>(sanisizer::attest_gez(data->rows())),
            sanisizer::cast<int>(sanisizer::attest_gez(data->cols())),
            static_cast<const double*>(data->begin())
        ));
    } else {
        auto mat = Rtatami::BoundNumericPointer(x);
        mptr.reset(new kmeans_tatami::Matrix<int, double, double, int>(mat->ptr));
    }

    const auto ndim = mptr->num_dimensions();
    const auto nobs = mptr->num_observations();
    auto centers = create_matrix<Rcpp::NumericMatrix>(ndim, nclusters);
    auto clusters = sanisizer::create<Rcpp::IntegerVector>(nobs);
    auto center_ptr = static_cast<double*>(centers.begin());
    auto cluster_ptr = static_cast<int*>(clusters.begin());

    std::unique_ptr<kmeans::Initialize<int, double, int, double> > iptr;
    if (init_method == "random") {
        auto ptr = new kmeans::InitializeRandom<int, double, int, double>;
        set_integer(random_seed, ptr->get_options().seed, "random.seed");
        iptr.reset(ptr);

    } else if (init_method == "kmeans++") {
        auto ptr = new kmeans::InitializeKmeanspp<int, double, int, double>;
        set_integer(nthreads, ptr->get_options().num_threads, "num.threads");
        set_integer(kmeanspp_seed, ptr->get_options().seed, "kmeanspp.seed"); 
        iptr.reset(ptr);

    } else if (init_method == "var-part") {
        auto ptr = new kmeans::InitializeVariancePartition<int, double, int, double>;
        set_bool(var_part_optimize_partition, ptr->get_options().optimize_partition, "var.part.optimize.partition");
        set_number(var_part_size_adjustment, ptr->get_options().size_adjustment, "var.part.size.adjustment");
        iptr.reset(ptr);

    } else {
        throw std::runtime_error("unknown init_method '" + init_method + "'");
    }

    std::unique_ptr<kmeans::Refine<int, double, int, double> > rptr;
    if (refine_method == "lloyd") {
        auto ptr = new kmeans::RefineLloyd<int, double, int, double>;
        set_integer(lloyd_iterations, ptr->get_options().max_iterations, "lloyd.iterations");
        set_integer(nthreads, ptr->get_options().num_threads, "num.threads"); 
        rptr.reset(ptr);

    } else if (refine_method == "hartigan-wong") {
        auto ptr = new kmeans::RefineHartiganWong<int, double, int, double>;
        set_integer(hartigan_wong_iterations, ptr->get_options().max_iterations, "hartigan.wong.iterations");
        set_integer(hartigan_wong_quick_transfer_iterations, ptr->get_options().max_quick_transfer_iterations, "hartigan.wong.quick.transfer.iterations");
        set_bool(hartigan_wong_quit_quick_transfer_failure, ptr->get_options().quit_on_quick_transfer_convergence_failure, "hartigan.wong.quit.quick.transfer.failure");
        set_integer(nthreads, ptr->get_options().num_threads, "num.threads");
        rptr.reset(ptr);
    }

    auto out = kmeans::compute(*mptr, *iptr, *rptr, nclusters, center_ptr, cluster_ptr);
    const auto actual_k = kmeans::remove_unused_centers(ndim, nobs, cluster_ptr, nclusters, center_ptr, out.sizes);
    if (actual_k != nclusters) {
        auto new_centers = create_matrix<Rcpp::NumericMatrix>(ndim, actual_k);
        std::copy_n(centers.begin(), new_centers.size(), new_centers.begin());
        centers = std::move(new_centers);
    }

    return Rcpp::List::create(
        Rcpp::Named("clusters") = clusters, 
        Rcpp::Named("centers") = centers,
        Rcpp::Named("iterations") = Rcpp::IntegerVector::create(out.iterations),
        Rcpp::Named("status") = Rcpp::IntegerVector::create(out.status) 
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_kmeans_defaults() {
    Rcpp::List output;

    {
        kmeans::InitializeRandomOptions opt;
        output["random.seed"] = opt.seed;
    }

    {
        kmeans::InitializeKmeansppOptions opt;
        output["num.threads"] = opt.num_threads; // just use this value as the representative, it should be 1 for everything anyway.
        output["kmeanspp.seed"] = opt.seed;
    }

    {
        kmeans::InitializeVariancePartitionOptions opt;
        output["var.part.optimize.partition"] = opt.optimize_partition;
        output["var.part.size.adjustment"] = opt.size_adjustment;
    }

    {
        kmeans::RefineLloydOptions opt;
        output["lloyd.iterations"] = opt.max_iterations;
    }

    {
        kmeans::RefineHartiganWongOptions opt;
        output["hartigan.wong.iterations"] = opt.max_iterations;
        output["hartigan.wong.quick.transfer.iterations"] = opt.max_quick_transfer_iterations;
        output["hartigan.wong.quit.quick.transfer.failure"] = opt.quit_on_quick_transfer_convergence_failure;
    }

    return output;
}
