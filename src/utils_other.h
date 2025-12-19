#ifndef UTILS_OTHER_H
#define UTILS_OTHER_H

#include "config.h"

#include <type_traits>

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

#endif
