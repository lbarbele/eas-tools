#ifndef _utl_vector_h
#define _utl_vector_h

#include <array>
#include <ostream>
#include <iomanip>
#include <cmath>

#include <util/type_traits.h>
#include <util/frame.h>

namespace util {

  class vector_base {};

  template <class T>
  class vector : vector_base {
  private:
    std::array<T, 3> m_data;
    frame_ptr m_frame;
  
  public:
    // - Constructors

    // * construct the zero vector
    vector(const frame_ptr frame = frame::standard): m_data{0, 0, 0}, m_frame(frame) {}

    // * construct vector with given coordinates
    vector(const T x, const T y, const T z, const frame_ptr frame = frame::standard) : m_data{x, y, z}, m_frame(frame) {}

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

    // - Vector normalization

    // * compute vector norm 
    auto norm() const
    {return std::sqrt((*this)*(*this));}

    // * normalize vector to given value and return it
    template <class U>
    auto& normalize(const U w = 1)
    {return (*this *= w/norm());}

    // - Frame-independent operations

    // * (assignment) multiplication by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector<T>& operator*=(const U scalar)
    {
      (*this)[0] *= scalar;
      (*this)[1] *= scalar;
      (*this)[2] *= scalar;
      return *this;
    }

    // * (assignment) division by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector<T>& operator/=(const U scalar)
    {
      (*this)[0] /= scalar;
      (*this)[1] /= scalar;
      (*this)[2] /= scalar;
      return *this;
    }

    // * multiplication by scalar
    template <class U, class R = decltype(T{} * U{}), typename = enable_if_scalar_t<U>>
    vector<R> operator*(const U scalar) const
    {
      return {
        (*this)[0] * scalar,
        (*this)[1] * scalar,
        (*this)[2] * scalar
      };
    }

    // * division by scalar
    template <class U, class R = decltype(T{} / U{}), typename = enable_if_scalar_t<U>>
    vector<R> operator/(const T scalar) const
    {
      return {
        (*this)[0] / scalar,
        (*this)[1] / scalar,
        (*this)[2] / scalar
      };
    }

    // * unary plus operator
    auto operator+() const
    {return (*this);}

    // * unary minus operator
    auto operator-() const
    {return (*this) * (-1.0);}

    // - Frame-dependent operations

    // * access the frame
    const frame_ptr& get_frame() const
    {return m_frame;}

    // * change frame
    vector<T>& set_frame(
      const frame_ptr& frame
    )
    {
      if (frame != get_frame()) {
        *this = frame->to() * get_frame()->from() * (*this);
        m_frame = frame;
      }
      return *this;
    }

    // * create a copy of vector in a different frame
    vector<T> on_frame(
      const frame_ptr& frame
    ) const
    {
      vector<T> other = (*this);
      other.set_frame(frame);
      return other;
    }

    // * (assignment) sum with vector
    template <class U, class R = decltype(T{} + U{})>
    vector<T>& operator+=(const vector<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      (*this)[0] += other[0];
      (*this)[1] += other[1];
      (*this)[2] += other[2];
      return (*this);
    }

    // * (assignment) subtraction with vector
    template <class U, class R = decltype(T{} - U{})>
    vector<T>& operator-=(const vector<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      (*this)[0] -= other[0];
      (*this)[1] -= other[1];
      (*this)[2] -= other[2];
      return (*this);
    }

    // * vector sum
    template <class U, class R = decltype(T{} + U{})>
    vector<T> operator+(const vector<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      return {
        (*this)[0] + other[0],
        (*this)[1] + other[1],
        (*this)[2] + other[2],
        get_frame()
      };
    }

    // * vector subtraction
    template <class U, class R = decltype(T{} + U{})>
    vector<T> operator-(const vector<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      return {
        (*this)[0] - other[0],
        (*this)[1] - other[1],
        (*this)[2] - other[2],
        get_frame()
      };
    }

    // * dot product
    template <class U, class R = decltype(T{} * U{})>
    R operator*(const vector<U>& v) const
    {
      const auto& other = v.on_frame(get_frame());
      return other[0]*(*this)[0] + other[1]*(*this)[1] + other[2]*(*this)[2];
    }

    template <class U, class R = decltype(T{} * U{})>
    vector<R> cross_product(
      const vector<U>& v
    ) const
    {
      other = v.on_frame(get_frame());
      return {
        (*this)[1]*other[2] - (*this)[2]*other[1],
        (*this)[2]*other[0] - (*this)[0]*other[2],
        (*this)[0]*other[1] - (*this)[1]*other[0]
      };
    }

  };

  // * scalar-vector product (from lhs)
  template <class T, class U, class R = decltype(T{} * U{}), typename = enable_if_scalar_t<T>>
  vector<R> operator*(
    const T scalar,
    const vector<U>& v
  )
  {
    return v*scalar;
  }

  // * matrix-vector product (from lhs)
  template <class T, class U, class R = decltype(T{} * U{})>
  vector<R> operator*(
    const square_matrix<T, 3>& mtx,
    const vector<U>& v
  )
  {
    return {
      mtx(0, 0)*v[0] + mtx(0, 1)*v[1] + mtx(0, 2)*v[2],
      mtx(1, 0)*v[0] + mtx(1, 1)*v[1] + mtx(1, 2)*v[2],
      mtx(2, 0)*v[0] + mtx(2, 1)*v[1] + mtx(2, 2)*v[2]
    };
  }

  // * print vector to output stream
  template<class T, class CharT, class Traits = std::char_traits<CharT> >
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const vector<T>& v
  )
  {
    auto w = stream.width() > 0? stream.width() : 15;
    for (auto & elem : v) {
      stream << std::setw(w) << elem;
    }
    return stream;
  }

  // - Aliases
  using vector_d = vector<double>;
  using vector_f = vector<float>;
}

#endif