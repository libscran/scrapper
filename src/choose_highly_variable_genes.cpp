//#include "config.h"

#include <vector>

#include "Rcpp.h"
#include "scran_variances/scran_variances.hpp"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector choose_highly_variable_genes(Rcpp::NumericVector stats, int top, bool larger, bool keep_ties, Rcpp::Nullable<Rcpp::NumericVector> bound) {
    scran_variances::ChooseHighlyVariableGenesOptions opt;
    opt.top = top;
    opt.larger = larger;
    opt.keep_ties = keep_ties;

    opt.use_bound = !bound.isNull();
    if (opt.use_bound) {
        Rcpp::NumericVector bnd(bound);
        if (bnd.size() != 1) {
            throw std::runtime_error("'bound' should be a numeric scalar");
        }
        opt.bound = bnd[0];
    }

    auto res = scran_variances::choose_highly_variable_genes_index(stats.size(), static_cast<const double*>(stats.begin()), opt);
    return Rcpp::IntegerVector(res.begin(), res.end());
}
