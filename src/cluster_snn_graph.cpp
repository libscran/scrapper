//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/scran_graph_cluster.hpp"

Rcpp::List cluster_multilevel(int num_vertices, Rcpp::IntegerVector edges, Rcpp::NumericVector weights, double resolution, int seed) {
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), static_cast<const int*>(edges.begin()), num_vertices, false);

    scran_graph_cluster::ClusterMultilevelOptions opt;
    opt.resolution = resolution;
    opt.seed = seed;

    scran_graph_cluster::ClusterMultilevelResults res;
    igraph_vector_t weight_view;
    igraph_vector_view(&weight_view, weights.begin(), weights.size());
    scran_graph_cluster::cluster_multilevel(graph.get(), &weight_view, opt, res);

    size_t nlevels = res.levels.ncol();
    Rcpp::List levels(nlevels);
    for (size_t l = 0; l < nlevels; ++l) {
        auto incol = res.levels.column(l);
        levels[l] = Rcpp::IntegerVector(incol.begin(), incol.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::IntegerVector::create(res.status),
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("levels") = levels,
        Rcpp::Named("modularity") = Rcpp::IntegerVector(res.modularity.begin(), res.modularity.end())
    );
}

Rcpp::List cluster_leiden(int num_vertices, Rcpp::IntegerVector edges, Rcpp::NumericVector weights, double resolution, bool use_cpm, int seed) {
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), static_cast<const int*>(edges.begin()), num_vertices, false);
    igraph_vector_t weight_view;
    igraph_vector_view(&weight_view, weights.begin(), weights.size());

    scran_graph_cluster::ClusterLeidenOptions opt;
    opt.resolution = resolution;
    opt.modularity = !use_cpm;
    opt.seed = seed;
    scran_graph_cluster::ClusterLeidenResults res;
    scran_graph_cluster::cluster_leiden(graph.get(), &weight_view, opt, res);

    return Rcpp::List::create(
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("quality") = Rcpp::NumericVector::create(res.quality)
    );
}

Rcpp::List cluster_leiden(int num_vertices, Rcpp::IntegerVector edges, Rcpp::NumericVector weights, int steps, int seed) {
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), static_cast<const int*>(edges.begin()), num_vertices, false);
    igraph_vector_t weight_view;
    igraph_vector_view(&weight_view, weights.begin(), weights.size());

    scran_graph_cluster::ClusterWalktrapOptions opt;
    opt.steps = steps;
    scran_graph_cluster::ClusterWalktrapResults res;
    scran_graph_cluster::cluster_walktrap(graph.get(), &weight_view, opt, res);

    size_t merge_ncol = res.merges.ncol();
    Rcpp::IntegerMatrix merges(res.merges.nrow(), merge_ncol);
    for (size_t m = 0; m < merge_ncol; ++m) {
        auto incol = res.merges.column(m);
        auto outcol = merges.column(m);
        std::copy(incol.begin(), incol.end(), outcol.begin());
    }

    return Rcpp::List::create(
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("merges") = merges,
        Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
    ); 
}
