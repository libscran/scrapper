#include "config.h"

#include <stdexcept>
#include <cstddef>

#include "scran_norm/scran_norm.hpp"
#include "sanisizer/sanisizer.hpp"

#include "utils_other.h"
#include "utils_block.h"

//[[Rcpp::export(rng=false)]]
Rcpp::List center_spike_in_factors(
    Rcpp::NumericVector endogenous,
    Rcpp::List spike_ins,
    Rcpp::Nullable<Rcpp::IntegerVector> block,
    bool lowest
) {
    auto block_info = MaybeBlock(block);
    auto ptr = block_info.get();
    const auto ncells = endogenous.size();
    auto e_output = Rcpp::clone(endogenous);

    const auto num_spikes = spike_ins.size();
    std::vector<Rcpp::NumericVector> s_output;
    sanisizer::reserve(s_output, num_spikes);
    std::vector<double*> s_ptrs;
    sanisizer::reserve(s_ptrs, num_spikes);

    for (I<decltype(num_spikes)> s = 0; s < num_spikes; ++s) {
        Rcpp::NumericVector current(spike_ins[s]);
        if (ncells != current.size()) {
            throw std::runtime_error("'endogenous' and 'spike.ins' should have the same number of cells");
        }
        s_output.push_back(Rcpp::clone(current));
        s_ptrs.push_back(static_cast<double*>(s_output.back().begin()));
    }

    if (ptr) {
        if (!sanisizer::is_equal(block_info.size(), ncells)) {
            throw std::runtime_error("'block' must be the same length as the number of cells");
        }

        scran_norm::CenterSpikeInFactorsBlockedOptions opt;
        opt.block_mode = (lowest ? scran_norm::CenterBlockMode::LOWEST : scran_norm::CenterBlockMode::PER_BLOCK);
        opt.ignore_invalid = true;
        scran_norm::center_spike_in_factors_blocked(
            sanisizer::cast<std::size_t>(ncells),
            static_cast<double*>(e_output.begin()),
            s_ptrs,
            ptr,
            opt
        );

    } else {
        scran_norm::CenterSpikeInFactorsOptions opt;
        opt.ignore_invalid = true;
        scran_norm::center_spike_in_factors(
            sanisizer::cast<std::size_t>(ncells),
            static_cast<double*>(e_output.begin()),
            s_ptrs,
            opt
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("endogenous") = e_output,
        Rcpp::Named("spike.ins") = Rcpp::List(s_output.begin(), s_output.end())
    );
}
