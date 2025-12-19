#ifndef UTILS_OTHER_H
#define UTILS_OTHER_H

#include <type_traits>

template<typename Input_>
using I = std::remove_reference_t<std::remove_cv_t<Input_> >;

#endif
