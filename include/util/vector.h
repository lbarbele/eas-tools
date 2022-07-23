#ifndef _util_vector_h
#define _util_vector_h

#include <array>
#include <ostream>
#include <iomanip>
#include <cmath>

#include <util/coordinates.h>
#include <util/type_traits.h>
#include <util/frame.h>

namespace util {

  template <class T>
  class vector_t : public vector_base_t, public coordinates_t<T> {
  private:
    frame_ptr m_frame;
  
  public:
    // - Constructors

    // * undefined coordinates on given frame (defaults to standard frame)
    vector_t(
      const frame_ptr frame = frame::standard
    ) :
      m_frame(frame)
    {}

    // * construct a vector from a coordinates object and frame (explicit)
    vector_t(
      const coordinates_t<T>& coordinates,
      const frame_ptr frame
    ) :
      coordinates_t<T>(coordinates),
      m_frame(frame)
    {}

    // * construct a point with explicit coordinates and frame (defaults to standard frame)
    vector_t(
      const T& x,
      const T& y,
      const T& z,
      const frame_ptr frame = frame::standard
    ) :
      vector_t({x, y, z}, frame)
    {}

    // - Vector normalization

    // * compute vector norm 
    auto norm() const
    {return std::hypot((*this)[0], (*this)[1], (*this)[2]);}

    // * normalize vector to given value and return it
    template <class U>
    auto& normalize(const U w = 1)
    {
      (*this) *= w/norm();
      return (*this);
    }

    // - Frame manipulation

    // * access the frame
    const frame_ptr& get_frame() const
    {return m_frame;}

    // * change frame
    vector_t<T>& set_frame(
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
    vector_t<T> on_frame(
      const frame_ptr& frame
    ) const
    {
      vector_t<T> other = (*this);
      other.set_frame(frame);
      return other;
    }

    // - Frame-independent operations

    // * (assignment) multiplication by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector_t<T>& operator*=(const U scalar)
    {
      (*this)[0] *= scalar;
      (*this)[1] *= scalar;
      (*this)[2] *= scalar;
      return *this;
    }

    // * (assignment) division by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector_t<T>& operator/=(const U scalar)
    {
      (*this)[0] /= scalar;
      (*this)[1] /= scalar;
      (*this)[2] /= scalar;
      return *this;
    }

    // * multiplication by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector_t<decltype(T{}*U{})> operator*(const U scalar) const
    {
      return {
        (*this)[0] * scalar,
        (*this)[1] * scalar,
        (*this)[2] * scalar,
        get_frame()
      };
    }

    // * division by scalar
    template <class U, typename = enable_if_scalar_t<U>>
    vector_t<decltype(T{}*U{})> operator/(const U scalar) const
    {
      return {
        (*this)[0] / scalar,
        (*this)[1] / scalar,
        (*this)[2] / scalar,
        get_frame()
      };
    }

    // * unary plus operator
    vector_t<decltype(+T{})> operator+() const
    {
      return {
        +(*this)[0],
        +(*this)[1],
        +(*this)[2],
        get_frame()
      };
    }

    // * unary minus operator
    vector_t<decltype(-T{})> operator-() const
    {
      return {
        -(*this)[0],
        -(*this)[1],
        -(*this)[2],
        get_frame()
      };
    }
    // - Frame-dependent operations

    // * (assignment) sum with vector
    template <class U>
    vector_t<T>& operator+=(const vector_t<U>& v)
    {
      const auto other = v.on_frame(get_frame());
      (*this)[0] += other[0];
      (*this)[1] += other[1];
      (*this)[2] += other[2];
      return (*this);
    }

    // * (assignment) subtraction with vector
    template <class U>
    vector_t<T>& operator-=(const vector_t<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      (*this)[0] -= other[0];
      (*this)[1] -= other[1];
      (*this)[2] -= other[2];
      return (*this);
    }

    // * vector sum
    template <class U>
    vector_t<decltype(T{}+U{})> operator+(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return {
        (*this)[0] + other[0],
        (*this)[1] + other[1],
        (*this)[2] + other[2],
        get_frame()
      };
    }

    // * vector subtraction
    template <class U>
    vector_t<decltype(T{}-U{})> operator-(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return {
        (*this)[0] - other[0],
        (*this)[1] - other[1],
        (*this)[2] - other[2],
        get_frame()
      };
    }

    // * dot product
    template <class U>
    decltype(T{}*U{}) operator*(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return other[0]*(*this)[0] + other[1]*(*this)[1] + other[2]*(*this)[2];
    }

    // * cross product
    template <class U>
    vector_t<decltype(T{}*U{})> cross_product(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return {
        (*this)[1]*other[2] - (*this)[2]*other[1],
        (*this)[2]*other[0] - (*this)[0]*other[2],
        (*this)[0]*other[1] - (*this)[1]*other[0]
      };
    }

  };

  // * scalar-vector product (from lhs)
  template <class T, class U, typename = enable_if_scalar_t<T>>
  vector_t<decltype(T{}*U{})> operator*(
    const T scalar,
    const vector_t<U>& v
  )
  {
    // scalar-vector product is commutative
    return v*scalar;
  }

  // - Aliases
  using vector_d = vector_t<double>;
  using vector_f = vector_t<float>;
}

#endif