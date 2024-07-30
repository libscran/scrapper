//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/scran_graph_cluster.hpp"
#include "tatami/tatami.hpp"

#include "utils_graph.h"

//[[Rcpp::export(rng=false)]]
SEXP build_snn_graph(Rcpp::IntegerMatrix neighbors, std::string scheme, int num_threads, bool raw) {
    const int* nptr = neighbors.begin();
    size_t nrow = neighbors.rows();
    scran_graph_cluster::BuildSnnGraphResults<igraph_integer_t, igraph_real_t> output;

    scran_graph_cluster::BuildSnnGraphOptions opt;
    opt.num_threads = num_threads;
    if (scheme == "ranked") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::RANKED;
    } else if (scheme == "number") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::NUMBER;
    } else if (scheme == "jaccard") {
        opt.weighting_scheme = scran_graph_cluster::SnnWeightScheme::JACCARD;
    } else {
        throw std::runtime_error("unknown weighting scheme '" + scheme + "'");
    }

    size_t ncells = neighbors.cols();
    scran_graph_cluster::build_snn_graph(
        ncells,
        [&](int i) -> tatami::ArrayView<int> {
            return tatami::ArrayView<int>(nptr + nrow * static_cast<size_t>(i), nrow);
        },
        [](int i) -> int {
            return i - 1;
        },
        opt,
        output
    );

    if (raw) {
        return BuildSnnGraphPointer(new decltype(output)(std::move(output)), true);

    } else {
        size_t nedges = output.edges.size();
        Rcpp::IntegerVector edges(nedges);
        int* eptr = edges.begin();
        for (size_t e = 0; e < nedges; ++e) {
            eptr[e] = output.edges[e] + 1; // get to 1-based indexing.
        }

        return Rcpp::List::create(
            Rcpp::Named("vertices") = Rcpp::IntegerVector::create(ncells),
            Rcpp::Named("edges") = edges,
            Rcpp::Named("weights") = Rcpp::NumericVector(output.weights.begin(), output.weights.end())
        );
    }
}
