//#include "config.h"

#include <vector>
#include <string>
#include <stdexcept>

#include "Rcpp.h"

#include "scran_markers/scran_markers.hpp"

//[[Rcpp::export(rng=false)]]
SEXP summarize_effects(int num_genes, int num_groups, Rcpp::NumericVector effects, int num_threads) {
    size_t expected = num_groups;
    expected *= num_groups;
    expected *= num_genes;
    if (effects.size() != expected) {
        throw std::runtime_error("'effects' does not have the expected length");
    }

    std::vector<Rcpp::NumericVector> min, mean, median, max;
    min.reserve(num_groups);
    mean.reserve(num_groups);
    median.reserve(num_groups);
    max.reserve(num_groups);
    std::vector<Rcpp::IntegerVector> min_rank;
    min_rank.reserve(num_groups);

    std::vector<scran_markers::SummaryBuffers<double, int> > groupwise;
    groupwise.resize(num_groups);
    for (int g = 0; g < num_groups; ++g) {
        min.emplace_back(num_genes);
        groupwise[g].min = min.back().begin();
        mean.emplace_back(num_genes);
        groupwise[g].mean = mean.back().begin();
        median.emplace_back(num_genes);
        groupwise[g].median = median.back().begin();
        max.emplace_back(num_genes);
        groupwise[g].max = max.back().begin();
        min_rank.emplace_back(num_genes);
        groupwise[g].min_rank = min_rank.back().begin();
    }

    scran_markers::SummarizeEffectsOptions opt;
    opt.num_threads = num_threads;
    scran_markers::summarize_effects(num_genes, num_groups, static_cast<const double*>(effects.begin()), groupwise, opt);

    Rcpp::List output(num_groups);
    for (int g = 0; g < num_groups; ++g) {
        output[g] = Rcpp::List::create(
            Rcpp::Named("min") = min[g],
            Rcpp::Named("mean") = mean[g],
            Rcpp::Named("median") = median[g],
            Rcpp::Named("max") = max[g],
            Rcpp::Named("min.rank") = min_rank[g]
        );
    }
    return output;
}
