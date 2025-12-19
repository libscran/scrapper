#ifndef UTILS_GRAPH_H
#define UTILS_GRAPH_H

#include "config.h"

#include <cstddef>
#include <vector>

#include "scran_graph_cluster/scran_graph_cluster.hpp"

struct GraphComponents {
    std::size_t vertices = 0;
    std::vector<igraph_integer_t> edges;
    bool weighted = false;
    std::vector<igraph_real_t> weights;
};

typedef Rcpp::XPtr<GraphComponents> GraphComponentsPointer;

#endif
