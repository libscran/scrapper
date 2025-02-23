//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "scran_graph_cluster/build_snn_graph.hpp"
#include "tatami/tatami.hpp"

#include "utils_graph.h"

//[[Rcpp::export(rng=false)]]
SEXP build_snn_graph(Rcpp::IntegerMatrix neighbors, std::string scheme, int num_threads) {
    const int* nptr = neighbors.begin();
    size_t nrow = neighbors.rows();

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
    scran_graph_cluster::BuildSnnGraphResults<igraph_integer_t, igraph_real_t> buffers;
    scran_graph_cluster::build_snn_graph(
        ncells,
        [&](int i) -> tatami::ArrayView<int> {
            return tatami::ArrayView<int>(nptr + nrow * static_cast<size_t>(i), nrow);
        },
        [](int i) -> int {
            return i - 1;
        },
        opt,
        buffers
    );

    GraphComponentsPointer output(new GraphComponents);
    output->vertices = ncells;
    output->weighted = true;
    output->edges = std::move(buffers.edges);
    output->weights = std::move(buffers.weights);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List graph_to_list(SEXP ptr0) {
    GraphComponentsPointer ptr(ptr0);
    const auto& edges = ptr->edges;

    size_t nedges = edges.size();
    Rcpp::IntegerVector edges_p1(nedges);
    for (size_t e = 0; e < nedges; ++e) {
        edges_p1[e] = edges[e] + 1; // get to 1-based indexing.
    }

    SEXP weights_copy = R_NilValue;
    if (ptr->weighted) {
        const auto& weights = ptr->weights;
        weights_copy = Rcpp::NumericVector(weights.begin(), weights.end());
    }

    return Rcpp::List::create(
        Rcpp::Named("vertices") = Rcpp::IntegerVector::create(ptr->vertices),
        Rcpp::Named("edges") = edges_p1,
        Rcpp::Named("weights") = weights_copy
    );
}

//[[Rcpp::export(rng=false)]]
SEXP list_to_graph(Rcpp::List contents) {
    GraphComponentsPointer output(new GraphComponents);

    if (contents.size() != 3) {
        throw std::runtime_error("'x' should be a list of length 3");
    }

    Rcpp::IntegerVector vertices(contents[0]);
    if (vertices.size() != 1 || vertices[0] < 0) {
        throw std::runtime_error("first element of 'x' should be a non-negative integer scalar");
    }
    output->vertices = vertices[0];

    Rcpp::IntegerVector edges_p1(contents[1]);
    size_t nedges = edges_p1.size();
    auto& edges = output->edges;
    edges.resize(nedges);
    for (size_t i = 0; i < nedges; ++i) {
        edges[i] = edges_p1[i] - 1; // get back to 0-based indexing.
    }

    Rcpp::RObject weights(contents[2]);
    if (!weights.isNULL()) {
        output->weighted = true;
        Rcpp::NumericVector weights_vec(weights);
        output->weights.insert(output->weights.end(), weights_vec.begin(), weights_vec.end());
    }

    return output;
}
