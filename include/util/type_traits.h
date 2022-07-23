#ifndef _util_type_traits_h
#define _util_type_traits_h

#include <type_traits>

namespace util {

  // * base classes
  class matrix_base_t {};
  class point_base_t {};
  class vector_base_t {};
  class coordinates_base_t {};

  // * enable if derivate from coordinates_base_t
  template <class T, class U = void>
  using enable_if_coordinates_t = std::enable_if_t<std::is_base_of_v<coordinates_base_t, T>, U>;

  // * enable if derived from point_base_t
  template <class T, class U = void>
  using enable_if_point_t = std::enable_if_t<std::is_base_of_v<point_base_t, T>, U>;

  // * enable if derived from vector_base_t
  template <class T, class U = void>
  using enable_if_vector_t = std::enable_if_t<std::is_base_of_v<vector_base_t, T>, U>;

  // * enable if derived from matrix_base_t
  template <class T, class U = void>
  using enable_if_matrix_t = std::enable_if_t<std::is_base_of_v<matrix_base_t, T>, U>;

  // * enable if not derived from coordinates_base_t nor from matrix
  template <class T, class U = void>
  using enable_if_scalar_t = std::enable_if_t<
    !std::is_base_of_v<matrix_base_t, T> && !std::is_base_of_v<coordinates_base_t, T>
  , U>;

}

#endif