//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/scran_graph_cluster.hpp"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::List build_snn_graph(Rcpp::IntegerMatrix neighbors, std::string scheme, int num_threads) {
    const int* nptr = neighbors.begin();
    size_t nrow = neighbors.rows();
    scran_graph_cluster::BuildSnnGraphResults<int, double> output;
    
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

    scran_graph_cluster::build_snn_graph(
        neighbors.cols(), 
        [&](int i) -> tatami::ArrayView<int> {
            return tatami::ArrayView<int>(nptr + nrow * static_cast<size_t>(i), nrow);
        },
        [](int i) -> int {
            return i;
        },
        opt,
        output
    );

    return Rcpp::List::create(
        Rcpp::Named("edges") = Rcpp::IntegerVector(output.edges.begin(), output.edges.end()),
        Rcpp::Named("weights") = Rcpp::NumericVector(output.weights.begin(), output.weights.end())
    );
}
