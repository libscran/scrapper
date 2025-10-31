//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/scran_graph_cluster.hpp"

#include "utils_graph.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_multilevel(SEXP ptr0, double resolution, int seed) {
    GraphComponentsPointer ptr(ptr0);
    const auto& edges = ptr->edges;

    raiigraph::initialize();
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), edges.data(), ptr->vertices, false); 

    igraph_vector_t* weight_view_ptr = NULL;
    igraph_vector_t weight_view{};
    if (ptr->weighted) {
        weight_view_ptr = &weight_view;
        weight_view = igraph_vector_view(ptr->weights.data(), ptr->weights.size());
    }

    scran_graph_cluster::ClusterMultilevelOptions opt;
    opt.resolution = resolution;
    opt.seed = seed;

    scran_graph_cluster::ClusterMultilevelResults res;
    scran_graph_cluster::cluster_multilevel(graph.get(), weight_view_ptr, opt, res);

    size_t nlevels = res.levels.nrow();
    Rcpp::List levels(nlevels);
    for (size_t l = 0; l < nlevels; ++l) {
        auto incol = res.levels.row(l);
        levels[l] = Rcpp::IntegerVector(incol.begin(), incol.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("levels") = levels,
        Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_leiden(SEXP ptr0, double resolution, std::string objective, int seed) {
    GraphComponentsPointer ptr(ptr0);
    const auto& edges = ptr->edges;

    raiigraph::initialize();
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), edges.data(), ptr->vertices, false); 

    igraph_vector_t weight_view{};
    igraph_vector_t* weight_view_ptr = NULL;
    if (ptr->weighted) {
        weight_view = igraph_vector_view(ptr->weights.data(), ptr->weights.size());
        weight_view_ptr = &weight_view;
    }

    scran_graph_cluster::ClusterLeidenOptions opt;
    opt.resolution = resolution;
    opt.seed = seed;
    opt.report_quality = true;

    if (objective == "modularity") {
        opt.objective = IGRAPH_LEIDEN_OBJECTIVE_MODULARITY;
    } else if (objective == "cpm") {
        opt.objective = IGRAPH_LEIDEN_OBJECTIVE_CPM;
    } else if (objective == "er") {
        opt.objective = IGRAPH_LEIDEN_OBJECTIVE_ER;
    } else {
        throw std::runtime_error("unknown Leiden objective '" + objective + "'");
    }

    scran_graph_cluster::ClusterLeidenResults res;
    scran_graph_cluster::cluster_leiden(graph.get(), weight_view_ptr, opt, res);

    return Rcpp::List::create(
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("quality") = Rcpp::NumericVector::create(res.quality)
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_walktrap(SEXP ptr0, int steps) {
    GraphComponentsPointer ptr(ptr0);
    const auto& edges = ptr->edges;

    raiigraph::initialize();
    auto graph = scran_graph_cluster::edges_to_graph(edges.size(), edges.data(), ptr->vertices, false); 

    igraph_vector_t* weight_view_ptr = NULL;
    igraph_vector_t weight_view{};
    if (ptr->weighted) {
        weight_view = igraph_vector_view(ptr->weights.data(), ptr->weights.size());
        weight_view_ptr = &weight_view;
    }

    scran_graph_cluster::ClusterWalktrapOptions opt;
    opt.steps = steps;

    scran_graph_cluster::ClusterWalktrapResults res;
    scran_graph_cluster::cluster_walktrap(graph.get(), weight_view_ptr, opt, res);

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
