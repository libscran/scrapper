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

inline int parse_single_integer(const Rcpp::IntegerVector& src, const std::string& arg) {
    if (src.size() != 1) {
        throw std::runtime_error("expected a single integer for '" + arg + "'");
    }
    return src[0];
}

inline bool parse_single_bool(const Rcpp::LogicalVector& src, const std::string& arg) {
    if (src.size() != 1) {
        throw std::runtime_error("expected a single boolean for '" + arg + "'");
    }
    if (src[0] == NA_LOGICAL) {
        throw std::runtime_error("expected a non-NA boolean for '" + arg + "'");
    }
    return src[0];
}

inline double parse_single_number(const Rcpp::NumericVector& src, const std::string& arg) {
    if (src.size() != 1) {
        throw std::runtime_error("expected a single number for '" + arg + "'");
    }
    return src[0];
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

inline std::string parse_single_string(const Rcpp::RObject& src, const std::string& arg) {
    if (src.sexp_type() != STRSXP) {
        throw std::runtime_error("expected a string for '" + arg + "'");
    }
    return parse_single_string(Rcpp::CharacterVector(src), arg);
}

template<typename Type_>
void set_integer(const Rcpp::RObject& source, Type_& target, const std::string& arg) {
    if (source.isNULL()) {
        return;
    }
    static_assert(std::is_integral<Type_>::value);
    if (source.sexp_type() == INTSXP) {
        target = sanisizer::cast<Type_>(parse_single_integer(Rcpp::IntegerVector(source), arg));
    } else if (source.sexp_type() == REALSXP) {
        target = sanisizer::from_float<Type_>(parse_single_number(Rcpp::NumericVector(source), arg));
    } else {
        throw std::runtime_error("expected an integer or NULL for '" + arg + "'");
    }
}

template<typename Type_>
void set_optional_integer(const Rcpp::RObject& source, std::optional<Type_>& target, const std::string& arg) {
    if (source.isNULL()) {
        return;
    }
    static_assert(std::is_integral<Type_>::value);
    if (source.sexp_type() == INTSXP) {
        target = sanisizer::cast<Type_>(parse_single_integer(Rcpp::IntegerVector(source), arg));
    } else if (source.sexp_type() == REALSXP) {
        target = sanisizer::from_float<Type_>(parse_single_number(Rcpp::NumericVector(source), arg));
    } else {
        throw std::runtime_error("expected an integer or NULL for '" + arg + "'");
    }
}

template<typename Type_>
void set_number(const Rcpp::RObject& source, Type_& target, const std::string& arg) {
    if (source.isNULL()) {
        return;
    }
    static_assert(std::is_floating_point<Type_>::value);
    if (source.sexp_type() == INTSXP) {
        target = parse_single_integer(Rcpp::IntegerVector(source), arg);
    } else if (source.sexp_type() == REALSXP) {
        target = parse_single_number(Rcpp::NumericVector(source), arg);
    } else {
        throw std::runtime_error("expected an integer or NULL for '" + arg + "'");
    }
}

template<typename Type_>
void set_optional_number(const Rcpp::RObject& source, std::optional<Type_>& target, const std::string& arg) {
    if (source.isNULL()) {
        return;
    }
    static_assert(std::is_floating_point<Type_>::value);
    if (source.sexp_type() == INTSXP) {
        target = parse_single_integer(Rcpp::IntegerVector(source), arg);
    } else if (source.sexp_type() == REALSXP) {
        target = parse_single_number(Rcpp::NumericVector(source), arg);
    } else {
        throw std::runtime_error("expected a number or NULL for '" + arg + "'");
    }
}

inline void set_bool(const Rcpp::RObject& source, bool& target, const std::string& arg) {
    if (source.isNULL()) {
        return;
    }
    if (source.sexp_type() == INTSXP) {
        target = parse_single_integer(Rcpp::IntegerVector(source), arg);
    } else if (source.sexp_type() == REALSXP) {
        target = parse_single_number(Rcpp::NumericVector(source), arg);
    } else if (source.sexp_type() == LGLSXP) {
        target = parse_single_bool(Rcpp::LogicalVector(source), arg);
    } else {
        throw std::runtime_error("expected a boolean or NULL for '" + arg + "'");
    }
}

#endif
