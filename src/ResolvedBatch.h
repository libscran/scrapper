#ifndef RESOLVED_BATCH_H
#define RESOLVED_BATCH_H

#include "Rcpp.h"
#include <algorithm>

class ResolvedBatch {
public:
    ResolvedBatch(Rcpp::Nullable<Rcpp::IntegerVector> batch) {
        my_has_batch = batch.isNotNull();
        if (my_has_batch) {
            my_batch = Rcpp::IntegerVector(batch);
        }
        return;
    }

    const int* get() const {
        return (my_has_batch ? static_cast<const int*>(my_batch.begin()) : NULL);
    }

    size_t number() const {
        size_t nbatches = 1;
        if (my_has_batch && my_batch.size()){
            nbatches = *std::max_element(my_batch.begin(), my_batch.end()) + 1;
        }
        return nbatches;
    }

private:
    bool my_has_batch;
    // Need to carry along the vector to avoid garbage 
    // collection of any allocated memory (e.g., ALTREP'd)
    // from the Rcpp::Integer initialization.
    Rcpp::IntegerVector my_batch;
};

#endif
