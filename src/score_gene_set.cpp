//#include "config.h"

#include <vector>
#include <stdexcept>

#include "Rcpp.h"
#include "gsdecon/gsdecon.hpp"
#include "Rtatami.h"

#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List score_gene_set(
    SEXP x,
    int rank,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    bool scale,
    bool realized,
    int irlba_work,
    int irlba_iterations,
    int irlba_seed,
    int num_threads)
{
    auto mat = Rtatami::BoundNumericPointer(x);
    const auto& matrix = *(mat->ptr);
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    gsdecon::Options opt;
    opt.rank = rank;
    opt.scale = scale;
    opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
    opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);
    opt.realize_matrix = realized;
    opt.irlba_options.extra_work = irlba_work;
    opt.irlba_options.max_iterations = irlba_iterations;
    opt.irlba_options.seed = irlba_seed;
    opt.num_threads = num_threads;

    size_t NR = matrix.nrow();
    size_t NC = matrix.ncol();
    Rcpp::NumericVector scores(NC), weights(NR);
    gsdecon::Buffers<double> output;
    output.scores = scores.begin();
    output.weights = weights.begin();

    if (ptr) {
        if (block_info.size() != NC) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }
        gsdecon::compute_blocked(matrix, block_info.get(), opt, output);

    } else {
        gsdecon::compute(matrix, opt, output);
    }

    return Rcpp::List::create(
        Rcpp::Named("scores") = scores,
        Rcpp::Named("weights") = weights
    );
}
