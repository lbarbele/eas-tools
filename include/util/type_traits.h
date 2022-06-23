#ifndef _util_type_traits_h
#define _util_type_traits_h

#include <type_traits>

namespace util {

  // * base classes
  class matrix_base;
  class vector_base;

  // * enable if derived from matrix_base
  template <class T, class U = void>
  using enable_if_matrix_t = std::enable_if_t<std::is_base_of_v<matrix_base, T>, U>;

  // * enable if derived from vector_base
  template <class T, class U = void>
  using enable_if_vector_t = std::enable_if_t<std::is_base_of_v<vector_base, T>, U>;

  // * enable if not derived from vector nor from matrix
  template <class T, class U = void>
  using enable_if_scalar_t = std::enable_if_t<
    !std::is_base_of_v<matrix_base, T> && !std::is_base_of_v<vector_base, T>
  , U>;

}

#endif