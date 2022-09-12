#ifndef _util_coordinates_h
#define _util_coordinates_h

#include <array>
#include <cmath>

#include <util/math.h>
#include <util/matrix.h>
#include <util/type_traits.h>
#include <util/units.h>

namespace util {

  // - implementation of a type holding raw coordinates

  template <concepts::scalar T>
  class coordinates_t : coordinates_base_t {
  public:
    using value_type = std::remove_cvref_t<T>;

  protected:
    std::array<value_type, 3> m_data;

  public:

    // * constructors

    coordinates_t() : m_data{value_type(0), value_type(0), value_type(0)} {}
    coordinates_t(const value_type& a, const value_type& b, const value_type& c) : m_data{a, b, c} {}

    template <std::convertible_to<value_type> U>
    coordinates_t(const coordinates_t<U> c) : coordinates_t(c.x(), c.y(), c.z()) {}

    // * compute spherical coordinates

    auto get_r() const
    {return util::math::hypot(x(), y(), z());}

    auto get_theta() const
    {
      const auto r = get_r();
      return r > value_type(1e-20)?
        util::math::acos(z()/r) :
        units::radian_t(0);
    }

    auto get_phi() const
    {
      return util::math::hypot(x(), y()) > value_type(1e-20)?
        util::math::atan2(y(), x()) :
        units::radian_t(0);
    }

    // * element access

    // named access to x, y, and z coordinates

    constexpr auto& x() {return m_data[0];}
    constexpr auto& y() {return m_data[1];}
    constexpr auto& z() {return m_data[2];}

    constexpr const auto& x() const {return m_data[0];}
    constexpr const auto& y() const {return m_data[1];}
    constexpr const auto& z() const {return m_data[2];}

    // indexed access with bounds checking

    constexpr auto& at(const size_t i)
    {return m_data.at(i);}

    constexpr const auto& at(const size_t i) const
    {return m_data.at(i);}

    // indexed access w/o bounds checking

    constexpr auto& operator[](const size_t i)
    {return m_data[i];}

    constexpr const auto& operator[](const size_t i) const
    {return m_data[i];}

    // * iterators

    constexpr auto begin() {return m_data.begin();}
    constexpr auto end() {return m_data.end();}
    
    constexpr auto begin() const {return m_data.begin();}
    constexpr auto cbegin() const {return m_data.cbegin();}
    constexpr auto end() const {return m_data.end();}
    constexpr auto cend() const {return m_data.cend();}

    constexpr auto rbegin() {return m_data.rbegin();}
    constexpr auto rend() {return m_data.rend();}

    constexpr auto rbegin() const {return m_data.rbegin();}
    constexpr auto crbegin() const {return m_data.crbegin();}
    constexpr auto rend() const {return m_data.rend();}
    constexpr auto crend() const {return m_data.crend();}
  };

  // * matrix transformation of coordinates for any type derived from coordinates_t

  template <class T, concepts::arithmetic U>
  requires std::is_base_of_v<coordinates_base_t, T>
  constexpr auto operator*(
    const square_matrix<U, 3>& mtx,
    const T& obj
  )
  {
    // copy every other data member of the derived type
    auto transformed = obj;

    // then transform the coordinates
    transformed.x() = mtx(0, 0)*obj.x() + mtx(0, 1)*obj.y() + mtx(0, 2)*obj.z();
    transformed.y() = mtx(1, 0)*obj.x() + mtx(1, 1)*obj.y() + mtx(1, 2)*obj.z();
    transformed.z() = mtx(2, 0)*obj.x() + mtx(2, 1)*obj.y() + mtx(2, 2)*obj.z();

    // done
    return transformed;
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
    stream << std::setw(0) << "(";
    stream << std::setw(w) << c.x() << ", ";
    stream << std::setw(w) << c.y() << ", ";
    stream << std::setw(w) << c.z() << ")";
    return stream;
  }


}

#endif