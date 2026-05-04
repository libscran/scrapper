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

template<typename Type_>
void set_integer(const Rcpp::Nullable<Rcpp::IntegerVector>& source, Type_& target, const std::string& arg) {
    if (source.isNull()) {
        return;
    }
    Rcpp::IntegerVector src(source);
    if (src.size() != 1) {
        throw std::runtime_error("expected an integer or NULL for '" + arg + "'");
    }
    static_assert(std::is_integral<Type_>::value);
    target = sanisizer::cast<Type_>(src[0]);
}

template<typename Type_>
void set_optional_integer(const Rcpp::Nullable<Rcpp::IntegerVector>& source, std::optional<Type_>& target, const std::string& arg) {
    if (source.isNull()) {
        return;
    }
    Rcpp::IntegerVector src(source);
    if (src.size() != 1) {
        throw std::runtime_error("expected an integer or NULL for '" + arg + "'");
    }
    static_assert(std::is_integral<Type_>::value);
    target = sanisizer::cast<Type_>(src[0]);
}

template<typename Type_>
void set_number(const Rcpp::Nullable<Rcpp::NumericVector>& source, Type_& target, const std::string& arg) {
    if (source.isNull()) {
        return;
    }
    Rcpp::NumericVector src(source);
    if (src.size() != 1) {
        throw std::runtime_error("expected a number or NULL for '" + arg + "'");
    }
    static_assert(std::is_floating_point<Type_>::value);
    target = src[0];
}

inline void set_bool(const Rcpp::Nullable<Rcpp::LogicalVector>& source, bool& target, const std::string& arg) {
    if (source.isNull()) {
        return;
    }
    Rcpp::LogicalVector src(source);
    if (src.size() != 1) {
        throw std::runtime_error("expected a boolean or NULL for '" + arg + "'");
    }
    target = src[0];
}

inline std::string parse_single_string(const Rcpp::CharacterVector& src, const std::string& arg) {
    if (src.size() != 1) {
        throw std::runtime_error("expected a single string for '" + arg + "'");
    }
    Rcpp::String first = src[0];
    if (first == NA_STRING) {
        throw std::runtime_error("expected a non-NA string for '" + arg + "'");
    }
    return first.get_cstring();
}

#endif
