//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/scran_graph_cluster.hpp"

#include "utils_graph.h"

std::pair<raiigraph::Graph, igraph_vector_t> formulate_graph(SEXP graph, std::vector<igraph_real_t>& wbuffer) {
    if (TYPEOF(graph) == VECSXP) {
        Rcpp::List contents(graph);
        if (contents.size() != 3) {
            throw std::runtime_error("'x' should be a list of length 3");
        }

        Rcpp::IntegerVector vertices(contents[0]);
        if (vertices.size() != 1 || vertices[0] < 0) {
            throw std::runtime_error("first element of 'x' should be an integer scalar");
        }

        Rcpp::IntegerVector edges(contents[1]);
        const int* eptr = edges.begin();
        size_t nedges = edges.size();
        std::vector<igraph_integer_t> minus_1(nedges);
        for (size_t i = 0; i < nedges; ++i) {
            minus_1[i] = eptr[i] - 1; // get back to 0-based indexing.
        }

        igraph_vector_t weight_view{};
        Rcpp::NumericVector weights(contents[2]);
        if constexpr(std::is_same<double, igraph_real_t>::value) {
            igraph_vector_view(&weight_view, weights.begin(), weights.size());
        } else {
            wbuffer.insert(wbuffer.end(), weights.begin(), weights.end());
            igraph_vector_view(&weight_view, wbuffer.data(), wbuffer.size());
        }

        return std::make_pair(
            scran_graph_cluster::edges_to_graph(edges.size(), minus_1.data(), vertices[0], false),
            weight_view
        );

    } else if (TYPEOF(graph) == EXTPTRSXP) {
        BuildSnnGraphPointer ptr(graph);
        igraph_vector_t weight_view{};
        igraph_vector_view(&weight_view, ptr->weights.data(), ptr->weights.size());
        return std::make_pair(
            scran_graph_cluster::convert_to_graph(*ptr),
            weight_view
        );

    } else {
        throw std::runtime_error("unsupported graph representation");
        return std::make_pair(raiigraph::Graph{}, igraph_vector_t{});
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_multilevel(SEXP graph, double resolution, int seed) {
    scran_graph_cluster::ClusterMultilevelOptions opt;
    opt.resolution = resolution;
    opt.seed = seed;

    std::vector<igraph_real_t> wbuffer;
    auto gpair = formulate_graph(graph, wbuffer);
    scran_graph_cluster::ClusterMultilevelResults res;
    scran_graph_cluster::cluster_multilevel(gpair.first.get(), &(gpair.second), opt, res);

    size_t nlevels = res.levels.nrow();
    Rcpp::List levels(nlevels);
    for (size_t l = 0; l < nlevels; ++l) {
        auto incol = res.levels.row(l);
        levels[l] = Rcpp::IntegerVector(incol.begin(), incol.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::IntegerVector::create(res.status),
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("levels") = levels,
        Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_leiden(SEXP graph, double resolution, bool use_cpm, int seed) {
    scran_graph_cluster::ClusterLeidenOptions opt;
    opt.resolution = resolution;
    opt.modularity = !use_cpm;
    opt.seed = seed;
    opt.report_quality = true;

    std::vector<igraph_real_t> wbuffer;
    auto gpair = formulate_graph(graph, wbuffer);
    scran_graph_cluster::ClusterLeidenResults res;
    scran_graph_cluster::cluster_leiden(gpair.first.get(), &(gpair.second), opt, res);

    return Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::IntegerVector::create(res.status),
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("quality") = Rcpp::NumericVector::create(res.quality)
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::List cluster_walktrap(SEXP graph, int steps) {
    scran_graph_cluster::ClusterWalktrapOptions opt;
    opt.steps = steps;

    std::vector<igraph_real_t> wbuffer;
    auto gpair = formulate_graph(graph, wbuffer);
    scran_graph_cluster::ClusterWalktrapResults res;
    scran_graph_cluster::cluster_walktrap(gpair.first.get(), &(gpair.second), opt, res);

    size_t merge_ncol = res.merges.ncol();
    Rcpp::IntegerMatrix merges(res.merges.nrow(), merge_ncol);
    for (size_t m = 0; m < merge_ncol; ++m) {
        auto incol = res.merges.column(m);
        auto outcol = merges.column(m);
        std::copy(incol.begin(), incol.end(), outcol.begin());
    }

    return Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::IntegerVector::create(res.status),
        Rcpp::Named("membership") = Rcpp::IntegerVector(res.membership.begin(), res.membership.end()),
        Rcpp::Named("merges") = merges,
        Rcpp::Named("modularity") = Rcpp::NumericVector(res.modularity.begin(), res.modularity.end())
    ); 
}
