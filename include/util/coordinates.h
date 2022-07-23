#ifndef _util_coordinates_h
#define _util_coordinates_h

#include <array>
#include <cmath>

#include <util/matrix.h>
#include <util/type_traits.h>

namespace util {

  template <class T>
  class coordinates_t : public coordinates_base_t
  {
  protected:
    std::array<T, 3> m_data;

  public:

    // - Constructors
    
    coordinates_t() {};
    coordinates_t(const T& x, const T& y, const T& z) : m_data{x, y, z} {}

    // - Destructor

    virtual ~coordinates_t() {}

    // - Spherical coordinates

    T get_r() const
    {
      return std::hypot(x(), y(), z());
    }

    T get_theta() const
    {
      const T r = get_r();
      return r > 1e-20? std::acos(z()/get_r()) : 0;
    }

    T get_phi() const
    {
      const T r_sin_theta = std::hypot(x(), y());
      return r_sin_theta > 1e-20? std::atan2(y(), x()) : 0;
    }

    // - Element access

    // * direct access to x, y, and z coordinates
    T& x() {return m_data[0];}
    T& y() {return m_data[1];}
    T& z() {return m_data[2];}

    const T& x() const {return m_data[0];}
    const T& y() const {return m_data[1];}
    const T& z() const {return m_data[2];}

    // * access with bounds checking
    T& at(const size_t i)
    {return m_data.at(i);}

    const T& at(const size_t i) const
    {return m_data.at(i);}

    // * no bounds checking
    T& operator[](const size_t i)
    {return m_data[i];}

    const T& operator[](const size_t i) const
    {return m_data[i];}

    // * standard iterators
    auto begin() {return m_data.begin();}
    auto end() {return m_data.end();}
    
    // * const iterators
    auto begin() const {return m_data.begin();}
    auto cbegin() const {return m_data.cbegin();}
    auto end() const {return m_data.end();}
    auto cend() const {return m_data.cend();}

    // * reverse iterators
    auto rbegin() {return m_data.rbegin();}
    auto rend() {return m_data.rend();}

    // * const reverse iterators
    auto rbegin() const {return m_data.rbegin();}
    auto crbegin() const {return m_data.crbegin();}
    auto rend() const {return m_data.rend();}
    auto crend() const {return m_data.crend();}
  };

  // * coordinates transformation for any (template) type derived from coordinates_t
  template <class T, class U, template<class> class CoordinatesT, typename = enable_if_coordinates_t<CoordinatesT<U>>>
  CoordinatesT<decltype(T{}*U{})>
  operator*(
    const square_matrix<T, 3>& mtx,
    const CoordinatesT<U>& c
  )
  {
    return {
      mtx(0, 0)*c.x() + mtx(0, 1)*c.y() + mtx(0, 2)*c.z(),
      mtx(1, 0)*c.x() + mtx(1, 1)*c.y() + mtx(1, 2)*c.z(),
      mtx(2, 0)*c.x() + mtx(2, 1)*c.y() + mtx(2, 2)*c.z()
    };
  }

  // * print coordinates to stream
  template<class T, class CharT, class Traits = std::char_traits<CharT> >
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const coordinates_t<T>& c
  )
  {
    auto w = stream.width() > 0? stream.width() : 15;
    stream << "(";
    stream << std::setw(w) << c.x() << ",";
    stream << std::setw(w) << c.y() << ",";
    stream << std::setw(w) << c.z() << ")";
    return stream;
  }


}

#endif