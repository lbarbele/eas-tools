#ifndef _util_vector_h
#define _util_vector_h

#include <array>
#include <ostream>
#include <iomanip>
#include <cmath>

#include <util/coordinates.h>
#include <util/type_traits.h>
#include <util/frame.h>
#include <util/math.h>

namespace util {

  // - forward declaration and aliases
  
  template <concepts::scalar T>
  class vector_t;

  using vector_d = vector_t<double>;
  using vector_f = vector_t<float>;

  // - implementation of the vector class

  template <concepts::scalar T>
  class vector_t : public coordinates_t<T> {
  public:
    using value_type = coordinates_t<T>::value_type;
    using frame_type = frame<value_type>;
    using frame_ptr_type = frame_type::ptr_type;

  private:
    frame_ptr_type m_frame;
  
  public:

    // - Constructors

    // * zero vector on given frame (defaults to standard frame)

    vector_t(
      const frame_ptr_type frm = frame_type::standard
    ) :
      coordinates_t<value_type>(value_type(0), value_type(0), value_type(0)),
      m_frame(frm)
    {}

    // * construct a vector from a coordinates object and frame (explicit)

    vector_t(
      const coordinates_t<value_type>& coordinates,
      const frame_ptr_type frm
    ) :
      coordinates_t<value_type>(coordinates),
      m_frame(frm)
    {}

    // * construct a point with explicit coordinates and frame (defaults to standard frame)

    vector_t(
      const value_type& x,
      const value_type& y,
      const value_type& z,
      const frame_ptr_type frm
    ) :
      vector_t({x, y, z}, frm)
    {}

    // * copy constructor from vector with compatible value type

    template <std::convertible_to<value_type> U>
    vector_t(
      const vector_t<U>& other
    ) :
      vector_t(other.x(), other.y(), other.z(), other.get_frame())
    {}

    // - Vector normalization

    // * compute vector norm 

    value_type norm() const
    {return this->get_r();}

    // * normalize vector to given value and return it

    vector_t& normalize(const value_type w)
    {
      (*this) *= w/norm();
      return (*this);
    }

    // * get copy of vector normalized to given value

    vector_t get_normalized(const value_type w) const
    {
      auto other = *this;
      other.normalize(w);
      return other;
    }

    // - Frame manipulation

    // * access the frame

    const frame_ptr_type& get_frame() const
    {return m_frame;}

    // * change frame

    vector_t& set_frame(
      const frame_ptr_type& frame
    )
    {
      if (frame != get_frame()) {
        *this = frame->to() * get_frame()->from() * (*this);
        m_frame = frame;
      }
      return *this;
    }

    // * create a copy of vector in a different frame

    vector_t on_frame(
      const frame_ptr_type& frame
    ) const
    {
      vector_t other = (*this);
      other.set_frame(frame);
      return other;
    }

    // * create a copy of the point on the frame of another object that has a frame

    auto on_frame_of(const auto& f) 
    -> vector_t<typename std::remove_cvref_t<decltype(f)>::value_type>
    {return on_frame(f.get_frame());}

    // - Frame-independent operations

    // * (assignment) multiplication by scalar

    vector_t& operator*=(const auto& x)
    requires std::convertible_to<decltype(value_type{}*x), value_type>
    {
      (*this)[0] *= x;
      (*this)[1] *= x;
      (*this)[2] *= x;
      return *this;
    }

    // * (assignment) division by scalar

    vector_t& operator/=(const auto& x)
    requires std::convertible_to<decltype(value_type{}/x), value_type>
    {
      (*this)[0] /= x;
      (*this)[1] /= x;
      (*this)[2] /= x;
      return *this;
    }

    // * multiplication by scalar

    auto operator*(const concepts::scalar auto& x) const
    {
      using ret_value_type = decltype(value_type{}*x);
      return vector_t<ret_value_type> {
        (*this)[0] * x,
        (*this)[1] * x,
        (*this)[2] * x,
        frame<ret_value_type>::create(get_frame()->to(), frame<ret_value_type>::standard)
      };
    }

    // * division by scalar

    auto operator/(const concepts::scalar auto& x) const
    {
      using ret_value_type = decltype(value_type{}/x);
      return vector_t<ret_value_type> {
        (*this)[0] / x,
        (*this)[1] / x,
        (*this)[2] / x,
        frame<ret_value_type>::create(get_frame()->to(), frame<ret_value_type>::standard)
      };
    }

    // * unary plus operator

    auto operator+() const
    {
      return vector_t<decltype(+value_type{})>{
        +(*this)[0],
        +(*this)[1],
        +(*this)[2],
        get_frame()
      };
    }

    // * unary minus operator

    auto operator-() const
    {
      return vector_t<decltype(-T{})> {
        -(*this)[0],
        -(*this)[1],
        -(*this)[2],
        get_frame()
      };
    }

    // - Frame-dependent operations

    // * (assignment) sum with vector

    template <std::convertible_to<value_type> U>
    vector_t& operator+=(const vector_t<U>& v)
    {
      const auto other = v.on_frame(get_frame());
      (*this)[0] += other[0];
      (*this)[1] += other[1];
      (*this)[2] += other[2];
      return (*this);
    }

    // * (assignment) subtraction with vector

    template <std::convertible_to<value_type> U>
    vector_t& operator-=(const vector_t<U>& v)
    {
      const auto& other = v.on_frame(get_frame());
      (*this)[0] -= other[0];
      (*this)[1] -= other[1];
      (*this)[2] -= other[2];
      return (*this);
    }

    // * vector sum

    template <std::common_with<value_type> U>
    auto operator+(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return vector_t<decltype(value_type{}+U{})> {
        (*this)[0] + other[0],
        (*this)[1] + other[1],
        (*this)[2] + other[2],
        get_frame()
      };
    }

    // * vector subtraction

    template <std::common_with<value_type> U>
    auto operator-(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return vector_t<decltype(value_type{}-U{})> {
        (*this)[0] - other[0],
        (*this)[1] - other[1],
        (*this)[2] - other[2],
        get_frame()
      };
    }

    // * dot product

    template <class U>
    auto operator*(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return other[0]*(*this)[0] + other[1]*(*this)[1] + other[2]*(*this)[2];
    }

    // * cross product

    template <class U>
    auto cross_product(
      const vector_t<U>& v
    ) const
    {
      using ret_value_type = decltype(value_type{}*U{});
      const auto other = v.on_frame(get_frame());
      return vector_t<ret_value_type> {
        (*this)[1]*other[2] - (*this)[2]*other[1],
        (*this)[2]*other[0] - (*this)[0]*other[2],
        (*this)[0]*other[1] - (*this)[1]*other[0],
        frame<ret_value_type>::create(get_frame()->to(), frame<ret_value_type>::standard)
      };
    }

    // * vector-vector comparison
    
    template <std::common_with<value_type> U>
    auto operator==(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return
        other.x() == this->x() &&
        other.y() == this->y() &&
        other.z() == this->z();
    }
  };

  // * scalar-vector product (from lhs)

  template <class T>
  auto operator*(const concepts::scalar auto& x, const vector_t<T> v)
  {return v * x;}

  // * compute angle or cos(angle) between two vectors

  template <class A, class B>
  auto cos_angle(const vector_t<A>& a, const vector_t<B>& b)
  {return (a * b) / (a.norm() * b.norm());}

  template <class A, class B>
  auto angle(const vector_t<A>& a, const vector_t<B>& b)
  {
    auto cosine = cos_angle(a, b);
    return math::acos(cosine);
  }
}

#endif