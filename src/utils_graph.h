#ifndef UTILS_GRAPH_H
#define UTILS_GRAPH_H

#include "Rcpp.h"

#include "scran_graph_cluster/scran_graph_cluster.hpp"

typedef Rcpp::XPtr<scran_graph_cluster::BuildSnnGraphResults<igraph_integer_t, igraph_real_t> > BuildSnnGraphPointer;

#endif
