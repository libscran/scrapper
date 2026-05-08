#include "config.h"

#include <stdexcept>
#include <string>

#include "gsdecon/gsdecon.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_other.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List score_gene_set(
    SEXP x,
    Rcpp::RObject rank,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject scale,
    Rcpp::RObject realized,
    Rcpp::RObject irlba_work,
    Rcpp::RObject irlba_iterations,
    Rcpp::RObject irlba_tolerance,
    Rcpp::RObject irlba_seed,
    Rcpp::RObject num_threads
) {
    auto mat = Rtatami::BoundNumericPointer(x);
    const auto& matrix = *(mat->ptr);
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();

    gsdecon::Options opt;
    set_integer(rank, opt.rank, "rank");
    set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
    set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");
    set_bool(scale, opt.scale, "scale");
    set_bool(realized, opt.realize_matrix, "realized");
    set_optional_integer(irlba_work, opt.irlba_options.extra_work, "extra.work");
    set_integer(irlba_iterations, opt.irlba_options.max_iterations, "iterations");
    set_number(irlba_tolerance, opt.irlba_options.convergence_tolerance, "tolerance");
    set_integer(irlba_seed, opt.irlba_options.seed, "seed");
    set_integer(num_threads, opt.num_threads, "num.threads");

    const auto NR = matrix.nrow();
    const auto NC = matrix.ncol();
    auto scores = sanisizer::create<Rcpp::NumericVector>(NC);
    auto weights = sanisizer::create<Rcpp::NumericVector>(NR);
    gsdecon::Buffers<double> output;
    output.scores = scores.begin();
    output.weights = weights.begin();

    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), NC)) {
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

//[[Rcpp::export(rng=false)]]
Rcpp::List score_gene_set_defaults() {
    Rcpp::List output;
    gsdecon::Options opt;
    output["rank"] = opt.rank;

    report_block_weight_policy_default(output, opt.block_weight_policy, "block.weight.policy", "scoreGeneSet");
    report_variable_block_weight_default(output, opt.variable_block_weight_parameters, "variable.block.weight");

    output["scale"] = opt.scale;
    output["realized"] = opt.realize_matrix;

    if (!opt.irlba_options.extra_work.has_value()) {
        output["extra.work"] = R_NilValue;
    } else {
        throw std::runtime_error("unexpected extra.work default for runPca"); 
    }

    output["iterations"] = opt.irlba_options.max_iterations;
    output["tolerance"] = opt.irlba_options.convergence_tolerance;
    output["seed"] = opt.irlba_options.seed;
    output["num.threads"] = opt.num_threads;
    return output;
}
