//#include "config.h"

#include "mumosa/mumosa.hpp"
#include "utils_block.h"
#include "BiocNeighbors.h"
#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector scale_by_neighbors(
    int num_cells,
    Rcpp::List embedding,
    int num_neighbors,
    Rcpp::Nullable<Rcpp::IntegerVector> block, 
    std::string block_weight_policy,
    Rcpp::NumericVector variable_block_weight,
    int num_threads,
    SEXP nn_builder
) {
    auto nmod = embedding.size();
    std::vector<std::pair<double, double> > values;
    values.reserve(nmod);

    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    BiocNeighbors::BuilderPointer builder(nn_builder);

    if (ptr) {
        mumosa::BlockedOptions opt;
        opt.num_neighbors = num_neighbors;
        opt.num_threads = num_threads;
        opt.block_weight_policy = parse_block_weight_policy(block_weight_policy);
        opt.variable_block_weight_parameters = parse_variable_block_weight(variable_block_weight);

        mumosa::BlockedIndicesFactory<int, int> factory(num_cells, ptr);
        auto buff = factory.create_buffers<double>();
        auto work = mumosa::create_workspace<double>(factory.sizes(), opt);

        std::vector<std::shared_ptr<const BiocNeighbors::Prebuilt> > prebuilts;
        for (decltype(nmod) x = 0; x < nmod; ++x) {
            Rcpp::NumericMatrix current(embedding[x]);
            factory.build(current.rows(), static_cast<const double*>(current.begin()), *builder, prebuilts, buff);
            values.push_back(mumosa::compute_distance_blocked(prebuilts, work, opt));
        }

    } else {
        auto dist = sanisizer::create<std::vector<double> >(num_cells); 
        mumosa::Options opt;
        opt.num_neighbors = num_neighbors;
        opt.num_threads = num_threads;

        for (decltype(nmod) x = 0; x < nmod; ++x) {
            Rcpp::NumericMatrix current(embedding[x]);
            const auto prebuilt = builder->build_unique(knncolle::SimpleMatrix(current.rows(), num_cells, static_cast<const double*>(current.begin())));
            values.push_back(mumosa::compute_distance<int, double>(*prebuilt, dist.data(), opt));
        }
    }

    auto output = mumosa::compute_scale<double>(values);
    return Rcpp::NumericVector(output.begin(), output.end());
}
