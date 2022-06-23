#ifndef _util_type_traits_h
#define _util_type_traits_h

#include <type_traits>

#include <util/matrix.h>
#include <util/vector.h>

namespace util {

  // * type T is derived from matrix_base
  template <class T, class U = void>
  using enable_if_matrix_t = std::enable_if_t<std::is_base_of_v<matrix_base, T>, U>;

  // * type T is derived from vector_base
  template <class T, class U = void>
  using enable_if_vector_t = std::enable_if_t<std::is_base_of_v<vector_base, T>, U>;

  // * type T is neither vector nor matrix
  template <class T, class U = void>
  using enable_if_scalar_t = std::enable_if_t<
    !std::is_base_of_v<matrix_base, T> && !std::is_base_of_v<vector_base, T>
  , U>;

}

#endif