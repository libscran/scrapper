#ifndef UTILS_GRAPH_H
#define UTILS_GRAPH_H

#include "Rcpp.h"

#include "scran_graph_cluster/scran_graph_cluster.hpp"

struct GraphComponents {
    size_t vertices;
    std::vector<igraph_integer_t> edges;
    bool weighted;
    std::vector<igraph_real_t> weights;
};

typedef Rcpp::XPtr<GraphComponents> GraphComponentsPointer;

#endif
