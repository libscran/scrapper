#ifndef UTILS_OTHER_H
#define UTILS_OTHER_H

#include "config.h"

#include <type_traits>
#include <optional>

#include "sanisizer/sanisizer.hpp"

template<typename Input_>
using I = std::remove_reference_t<std::remove_cv_t<Input_> >;

template<typename Matrix_, typename Rows_, typename Cols_>
Matrix_ create_matrix(Rows_ rows, Cols_ cols) {
    return Matrix_(
        sanisizer::cast<I<decltype(std::declval<Matrix_>().nrow())> >(rows),
        sanisizer::cast<I<decltype(std::declval<Matrix_>().ncol())> >(cols)
    );
}

template<typename Vector_, typename Type_>
auto set_optional_integer(const Rcpp::Nullable<Vector_>& source, std::optional<Type_>& target) {
    if (source.isNull()) {
        return;
    }

    static_assert(std::is_same<Vector_, Rcpp::IntegerVector>::value);
    static_assert(std::is_integral<Type_>::value);
    Rcpp::IntegerVector src(source);
    if (src.size() != 1) {
        throw std::runtime_error("expected an integer scalar or NULL");
    }
    target = sanisizer::cast<Type_>(src[0]);
}

#endif
