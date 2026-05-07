#include "config.h"

#include <vector>
#include <string>
#include <memory>

#include "mumosa/mumosa.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_block.h"
#include "utils_other.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector scale_by_neighbors(
    int num_cells,
    Rcpp::List embedding,
    Rcpp::RObject num_neighbors,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    Rcpp::RObject block_weight_policy,
    Rcpp::RObject variable_block_weight,
    Rcpp::RObject num_threads,
    SEXP nn_builder
) {
    const auto nmod = embedding.size();
    std::vector<std::pair<double, double> > values;
    sanisizer::reserve(values, nmod);

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    BiocNeighbors::BuilderPointer builder(nn_builder);

    if (ptr) {
        mumosa::BlockedOptions opt;
        set_integer(num_neighbors, opt.num_neighbors, "num.neighbors");
        set_integer(num_threads, opt.num_threads, "num.threads");
        set_block_weight_policy(block_weight_policy, opt.block_weight_policy, "block.weight.policy");
        set_variable_block_weight(variable_block_weight, opt.variable_block_weight_parameters, "variable.block.weight");

        mumosa::BlockedIndicesFactory<int, int> factory(num_cells, ptr);
        auto buff = factory.create_buffers<double>();
        auto work = mumosa::create_workspace<double>(factory.sizes(), opt);

        std::vector<std::shared_ptr<const BiocNeighbors::Prebuilt> > prebuilts;
        for (I<decltype(nmod)> x = 0; x < nmod; ++x) {
            Rcpp::NumericMatrix current(embedding[x]);
            factory.build(  
                sanisizer::cast<std::size_t>(current.rows()),
                static_cast<const double*>(current.begin()),
                *builder,
                prebuilts,
                buff
            );
            values.push_back(mumosa::compute_distance_blocked(prebuilts, work, opt));
        }

    } else {
        auto dist = sanisizer::create<std::vector<double> >(num_cells); 
        mumosa::Options opt;
        set_integer(num_neighbors, opt.num_neighbors, "num.neighbors");
        set_integer(num_threads, opt.num_threads, "num.threads");

        for (I<decltype(nmod)> x = 0; x < nmod; ++x) {
            Rcpp::NumericMatrix current(embedding[x]);
            const auto prebuilt = builder->build_unique(
                knncolle::SimpleMatrix(sanisizer::cast<std::size_t>(current.rows()),
                num_cells,
                static_cast<const double*>(current.begin()))
            );
            values.push_back(mumosa::compute_distance<int, double>(*prebuilt, dist.data(), opt));
        }
    }

    auto output = mumosa::compute_scale<double>(values);
    return Rcpp::NumericVector(output.begin(), output.end());
}

//[[Rcpp::export(rng=false)]]
Rcpp::List scale_by_neighbors_defaults(bool use_block) {
    Rcpp::List output;

    if (use_block) {
        mumosa::BlockedOptions opt;
        output["num.neighbors"] = opt.num_neighbors;
        output["num.threads"] = opt.num_threads;
    } else {
        mumosa::Options opt;
        output["num.neighbors"] = opt.num_neighbors;
        output["num.threads"] = opt.num_threads;
    }

    // Setting this for consistency, regardless of whether use_block = true or not..
    mumosa::BlockedOptions opt;
    report_block_weight_policy_default(output, opt.block_weight_policy, "block.weight.policy", "scaleByNeighbors");
    report_variable_block_weight_default(output, opt.variable_block_weight_parameters, "variable.block.weight");

    return output;
}
