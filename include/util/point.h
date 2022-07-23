#ifndef _util_point_h
#define _util_point_h

#include <array>
#include <cmath>

#include <util/frame.h>
#include <util/vector.h>

namespace util {

  class point_base_t {};

  template <class T>
  class point_t : public point_base_t {
  private:
    std::array<T, 3> m_data;
    frame_ptr m_frame;

  public:

    // - Constructors

    // * construct a point with 0 coordinates on given frame  
    point_t(const frame_ptr frame = frame::standard): m_data{0, 0, 0}, m_frame(frame) {}

    // * constructe a point with given coordinates and frame
    point_t(const T x, const T y, const T z, const frame_ptr frame = frame::standard) : m_data{x, y, z}, m_frame(frame) {}

    // - Element access

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

    // * direct access to x, y, and z coordinates
    T& x = m_data[0];
    T& y = m_data[1];
    T& z = m_data[2];

    // - Frame manipulation

    // * access the frame
    const frame_ptr& get_frame() const
    {return m_frame;}

    // * change the frame
    point_t<T>& set_frame(
      const frame_ptr& frame
    )
    {
      if (frame != get_frame()) {
        // rotate to the standard frame
        *this = get_frame()->from() * (*this);
        // move origin to the new frame
        x += get_frame()->origin().x - frame->origin().x;
        y += get_frame()->origin().y - frame->origin().y;
        z += get_frame()->origin().z - frame->origin().z;
        // rotate to the new frame
        *this = frame->to() * (*this);
        // set the frame
        m_frame = frame;
      }
      return *this;
    }

    // * create a copy of this point in a different frame
    point_t<T> on_frame(
      const frame_ptr& frame
    )
    {
      point_t other = (*this);
      other.set_frame(frame);
      return other;
    }

    // - Operations

    // * distance to origin and distance to point
    T distance() const
    {return std::hypot(x, y, z);}

    template <class U>
    decltype(T{} - U{}) distance(point_t<U> p) const
    {
      p.on_frame(get_frame());
      return std::hypot(x-p.x, y-p.y, z-p.z);
    }

    // * point-point subtraction produces a vector
    template <class U>
    vector<decltype(T{}-U{})>
    operator-(
      point_t<U> p
    ) const
    {
      p.set_frame(get_frame());
      return {x - p.x, y - p.y, z - p.z, get_frame()};
    }

    // * point-vector sum/subtraction produces a point
    template <class U>
    point_t<decltype(T{}+U{})>
    operator+(
      vector<U> v
    ) const
    {
      v.set_frame(get_frame());
      return {x + v[0], y + v[1], z + v[2], get_frame()};
    }

    template <class U>
    point_t<T>&
    operator+=(
      vector<U> v
    )
    {
      v.set_frame(get_frame());
      x += v[0];
      y += v[1];
      z += v[2];
      return *this;
    }

    template <class U>
    point_t<decltype(T{}-U{})>
    operator-(
      vector<U> v
    ) const
    {
      v.set_frame(get_frame());
      return {x - v[0], y - v[1], z - v[2], get_frame()};
    }

    template <class U>
    point_t<T>&
    operator-=(
      vector<U> v
    )
    {
      v.set_frame(get_frame());
      x -= v[0];
      y -= v[1];
      z -= v[2];
      return *this;
    }

  };

  // * product with a matrix (transformations) from lhs
  template <class T, class U, class R = decltype(T{} * U{})>
  point_t<R> operator*(
    const square_matrix<T, 3>& mtx,
    const point_t<U>& p
  )
  {
    return {
      mtx(0, 0)*p.x + mtx(0, 1)*p.y + mtx(0, 2)*p.z,
      mtx(1, 0)*p.x + mtx(1, 1)*p.y + mtx(1, 2)*p.z,
      mtx(2, 0)*p.x + mtx(2, 1)*p.y + mtx(2, 2)*p.z
    };
  }

  // * print point to output stream
  template<class T, class CharT, class Traits = std::char_traits<CharT> >
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const point_t<T>& p
  )
  {
    auto w = stream.width() > 0? stream.width() : 15;
    for (const auto& coord : p) {
      stream << std::setw(w) << coord;
    }
    return stream;
  }

  // - Aliases
  using point_d = point_t<double>;
  using point_f = point_t<float>;

}

#endif