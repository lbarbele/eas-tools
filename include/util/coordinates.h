#ifndef _util_coordinates_h
#define _util_coordinates_h

#include <array>
#include <cmath>

#include <util/math.h>
#include <util/matrix.h>
#include <util/type_traits.h>
#include <util/units.h>

namespace util {

  template <concepts::scalar T>
  struct coordinates_t : public coordinates_base_t {

    // - Type aliases
    
    using value_type = std::remove_cvref_t<T>;

    // - Data

    std::array<value_type, 3> m_data;

    // - Constructors
    
    constexpr coordinates_t() {};

    constexpr coordinates_t(
      const value_type& x,
      const value_type& y,
      const value_type& z
    ) :
      m_data{x, y, z}
    {}

    // - Destructor

    virtual ~coordinates_t() {}

    // - Spherical coordinates

    constexpr auto get_r() const
    {return util::math::hypot(x(), y(), z());}

    constexpr auto get_theta() const
    {
      const auto r = get_r();
      return r > value_type(1e-20)?
        util::math::acos(z()/r) :
        units::radian_t(0);
    }

    constexpr auto get_phi() const
    {
      return util::math::hypot(x(), y()).get_value() > 1e-20?
        util::math::atan2(y(), x()) :
        units::radian_t(0);
    }

    // - Element access

    // * direct access to x, y, and z coordinates

    constexpr auto& x() {return m_data[0];}
    constexpr auto& y() {return m_data[1];}
    constexpr auto& z() {return m_data[2];}

    constexpr const auto& x() const {return m_data[0];}
    constexpr const auto& y() const {return m_data[1];}
    constexpr const auto& z() const {return m_data[2];}

    // * access with bounds checking

    constexpr auto& at(const size_t i)
    {return m_data.at(i);}

    constexpr const auto& at(const size_t i) const
    {return m_data.at(i);}

    // * no bounds checking

    constexpr auto& operator[](const size_t i)
    {return m_data[i];}

    constexpr const auto& operator[](const size_t i) const
    {return m_data[i];}

    // * standard iterators

    constexpr auto begin() {return m_data.begin();}
    constexpr auto end() {return m_data.end();}
    
    // * const iterators

    constexpr auto begin() const {return m_data.begin();}
    constexpr auto cbegin() const {return m_data.cbegin();}
    constexpr auto end() const {return m_data.end();}
    constexpr auto cend() const {return m_data.cend();}

    // * reverse iterators

    constexpr auto rbegin() {return m_data.rbegin();}
    constexpr auto rend() {return m_data.rend();}

    // * const reverse iterators

    constexpr auto rbegin() const {return m_data.rbegin();}
    constexpr auto crbegin() const {return m_data.crbegin();}
    constexpr auto rend() const {return m_data.rend();}
    constexpr auto crend() const {return m_data.crend();}
  };

  // * coordinates transformation for any (template) type derived from coordinates_t

  template <class T, class U, template<class> class CoordT>
  requires std::derived_from<CoordT<U>, coordinates_t<U>>
  constexpr auto operator*(
    const square_matrix<T, 3>& mtx,
    const CoordT<U>& c
  )
  {
    CoordT<decltype(T{}*U{})> other = c;
    other.x() = mtx(0, 0)*c.x() + mtx(0, 1)*c.y() + mtx(0, 2)*c.z();
    other.y() = mtx(1, 0)*c.x() + mtx(1, 1)*c.y() + mtx(1, 2)*c.z();
    other.z() = mtx(2, 0)*c.x() + mtx(2, 1)*c.y() + mtx(2, 2)*c.z();
    return other;
  }

  // * print coordinates to stream
  template<class T, class CharT, class Traits>
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const coordinates_t<T>& c
  )
  {
    auto w = stream.width();
    stream << "(";
    stream << std::setw(w) << c.x() << ", ";
    stream << std::setw(w) << c.y() << ", ";
    stream << std::setw(w) << c.z() << ")";
    return stream;
  }


}

#endif