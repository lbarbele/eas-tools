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

  // - implementation of the vector_t class

  template <concepts::scalar T>
  class vector_t : public coordinates_t<T> {
  public:
    using value_type = coordinates_t<T>::value_type;

  private:
    frame_ptr m_frame;
  
  public:

    // * Constructors

    // zero vector on given frame (defaults to standard frame)

    vector_t(
      const frame_ptr frm = frame::standard
    ) :
      coordinates_t<value_type>{value_type(0), value_type(0), value_type(0)},
      m_frame(frm)
    {}

    // construct a vector from a coordinates object and frame

    vector_t(
      const coordinates_t<value_type>& coordinates,
      const frame_ptr frm = frame::standard
    ) :
      coordinates_t<value_type>(coordinates),
      m_frame(frm)
    {}

    // construct a vector with explicit coordinates and frame

    vector_t(
      const value_type& x,
      const value_type& y,
      const value_type& z,
      const frame_ptr frm = frame::standard
    ) :
      vector_t({x, y, z}, frm)
    {}

    // copy constructor from vector with compatible value type

    template <std::convertible_to<value_type> U>
    vector_t(
      const vector_t<U>& other
    ) :
      vector_t(other.x(), other.y(), other.z(), other.get_frame())
    {}

    // * Vector normalization

    // compute vector norm 

    value_type norm() const
    {return this->get_r();}

    // normalize vector to given value and return it

    vector_t& normalize(const value_type w)
    {return (*this) *= w/norm();}

    // get copy of vector normalized to given value

    auto get_normalized(const concepts::scalar auto w) const
    {return (*this / norm()) * w;}

    // * Frame manipulation

    // access the frame

    const auto& get_frame() const
    {return m_frame;}

    // change the frame

    vector_t& set_frame(const frame_ptr& frm)
    {
      if (frm != get_frame()) {
        *this = frm->to() * get_frame()->from() * (*this);
        m_frame = frm;
      }
      return *this;
    }

    // create a copy of this vector in a different frame

    vector_t on_frame(const frame_ptr& frame) const
    {
      vector_t other = (*this);
      return other.set_frame(frame);
    }

    // create a copy of the point on the frame of another object that has a frame

    vector_t on_frame_of(const auto& f) const
    {return on_frame(f.get_frame());}

    // * Frame-independent operations

    // (assignment) multiplication by scalar

    vector_t& operator*=(const concepts::scalar auto& x)
    requires std::convertible_to<decltype(value_type{}*x), value_type>
    {
      (*this)[0] *= x;
      (*this)[1] *= x;
      (*this)[2] *= x;
      return *this;
    }

    // (assignment) division by scalar

    vector_t& operator/=(const concepts::scalar auto& x)
    requires std::convertible_to<decltype(value_type{}/x), value_type>
    {
      (*this)[0] /= x;
      (*this)[1] /= x;
      (*this)[2] /= x;
      return *this;
    }

    // multiplication by scalar

    auto operator*(const concepts::scalar auto& x) const
    {
      return vector_t<decltype(value_type{}*x)> {
        (*this)[0] * x,
        (*this)[1] * x,
        (*this)[2] * x,
        get_frame()
      };
    }

    // division by scalar

    auto operator/(const concepts::scalar auto& x) const
    {
      return vector_t<decltype(value_type{}/x)> {
        (*this)[0] / x,
        (*this)[1] / x,
        (*this)[2] / x,
        get_frame()
      };
    }

    // unary plus operator

    auto operator+() const
    {
      return vector_t<decltype(+value_type{})>{
        +(*this)[0],
        +(*this)[1],
        +(*this)[2],
        get_frame()
      };
    }

    // unary minus operator

    auto operator-() const
    {
      return vector_t<decltype(-value_type{})> {
        -(*this)[0],
        -(*this)[1],
        -(*this)[2],
        get_frame()
      };
    }

    // * Frame-dependent operations

    // (assignment) sum with vector

    template <std::convertible_to<value_type> U>
    vector_t& operator+=(const vector_t<U>& v)
    {
      const auto other = v.on_frame_of(*this);
      (*this)[0] += other[0];
      (*this)[1] += other[1];
      (*this)[2] += other[2];
      return (*this);
    }

    // (assignment) subtraction with vector

    template <std::convertible_to<value_type> U>
    vector_t& operator-=(const vector_t<U>& v)
    {
      const auto& other = v.on_frame_of(*this);
      (*this)[0] -= other[0];
      (*this)[1] -= other[1];
      (*this)[2] -= other[2];
      return (*this);
    }

    // vector sum

    template <std::common_with<value_type> U>
    auto operator+(const vector_t<U>& v) const
    {
      const auto other = v.on_frame_of(*this);
      return vector_t<decltype(value_type{}+U{})> {
        (*this)[0] + other[0],
        (*this)[1] + other[1],
        (*this)[2] + other[2],
        get_frame()
      };
    }

    // vector subtraction

    template <std::common_with<value_type> U>
    auto operator-(const vector_t<U>& v) const
    {
      const auto other = v.on_frame_of(*this);
      return vector_t<decltype(value_type{}-U{})> {
        (*this)[0] - other[0],
        (*this)[1] - other[1],
        (*this)[2] - other[2],
        get_frame()
      };
    }

    // dot product

    template <concepts::scalar U>
    auto operator*(const vector_t<U>& v) const
    {
      const auto other = v.on_frame(get_frame());
      return other.x()*this->x() + other.y()*this->y() + other.z()*this->z();
    }

    // cross product

    template <concepts::scalar U>
    auto cross_product(const vector_t<U>& v) const
    {
      const auto other = v.on_frame_of(*this);
      return vector_t<decltype(value_type{}*U{})> {
        this->y()*other.z() - this->z()*other.y(),
        this->z()*other.x() - this->x()*other.z(),
        this->x()*other.y() - this->y()*other.x(),
        get_frame()
      };
    }

    // vector-vector comparison
    
    template <std::common_with<value_type> U>
    auto operator==(const vector_t<U>& v) const
    {return this->m_data == v.on_frame_of(*this).m_data;}
  };

  // * scalar-vector product (from lhs)

  template <concepts::scalar T>
  auto operator*(const concepts::scalar auto& x, const vector_t<T> v)
  {return v * x; } // commutative!

  // * compute angle or cos(angle) between two vectors

  template <concepts::scalar A, concepts::scalar B>
  auto cos_angle(const vector_t<A>& a, const vector_t<B>& b)
  {return a.get_normalized(1) * b.get_normalized(1);}

  template <concepts::scalar A, concepts::scalar B>
  auto angle(const vector_t<A>& a, const vector_t<B>& b)
  {return math::acos(cos_angle(a, b));}
}

#endif