#ifndef _util_point_h
#define _util_point_h

#include <array>
#include <cmath>
#include <compare>
#include <concepts>
#include <type_traits>

#include <util/coordinates.h>
#include <util/frame.h>
#include <util/vector.h>

namespace util {

  // - forward declarations and aliases

  template <concepts::scalar T>
  class point_t;

  using point_d = point_t<double>;
  using point_f = point_t<float>;

  // - implementation of the point_t class

  template <concepts::scalar T>
  class point_t : public point_base_t, public coordinates_t<T> {
  public:
    using value_type = coordinates_t<T>::value_type;
    using frame_type = frame<value_type>;
    using frame_ptr_type = frame_type::ptr_type;

  private:
    frame_ptr_type m_frame = frame_type::standard;

  public:

    // - Constructors

    // * origin of given frame (defaults to standard frame)

    point_t(
      const frame_ptr_type frm = frame_type::standard
    ) :
      coordinates_t<value_type>(value_type(0), value_type(0), value_type(0)),
      m_frame(frm)
    {}

    // * construct a point from a coordinates object and frame (explicit)

    point_t(
      const coordinates_t<value_type>& coordinates,
      const frame_ptr_type frm
    ) :
      coordinates_t<value_type>(coordinates),
      m_frame(frm)
    {}

    // * construct a point with given coordinates and frame

    point_t(
      const value_type& x,
      const value_type& y,
      const value_type& z,
      const frame_ptr_type frm = frame_type::standard
    ) :
      point_t({x, y, z}, frm)
    {}

    // * copy constructor from point with compatible value type

    template <std::convertible_to<value_type> U>
    point_t(
      const point_t<U>& other
    ) :
      point_t(other.x(), other.y(), other.z(), other.get_frame())
    {}


    // - Frame manipulation

    // * access the frame

    const auto& get_frame() const
    {return m_frame;}

    // * change the frame

    point_t& set_frame(
      const frame_ptr_type& frame
    )
    {
      if (frame != get_frame()) {
        // copy coordinates of this point but rotated to coincide with the std frame
        coordinates_t<value_type> new_coordinates = get_frame()->from() * (*this);
        // move origin of the coordinates to the new frame
        new_coordinates.x() += get_frame()->origin()[0] - frame->origin()[0];
        new_coordinates.y() += get_frame()->origin()[1] - frame->origin()[1];
        new_coordinates.z() += get_frame()->origin()[2] - frame->origin()[2];
        // rotate coordinates to the new frame
        new_coordinates = frame->to() * new_coordinates;
        // set new frame and coordinates of this point
        this->x() = new_coordinates.x();
        this->y() = new_coordinates.y();
        this->z() = new_coordinates.z();
        m_frame = frame;
      }
      return *this;
    }

    // * create a copy of this point in a different frame

    point_t on_frame(
      const frame_ptr_type& frame
    ) const
    {
      point_t other = *this;
      other.set_frame(frame);
      return other;
    }

    // * create a copy of the point on the frame of another object that has a frame

    auto on_frame_of(const auto& f) 
    -> point_t<typename std::remove_cvref_t<decltype(f)>::value_type>
    {return on_frame(f.get_frame());}

    // - Operations

    // * distance to another point

    template <std::common_with<value_type> U>
    value_type distance(const point_t<U>& p) const
    {return (p - *this).norm();}

    // * point-point subtraction produces a vector

    template <std::common_with<value_type> U>
    auto operator-(
      const point_t<U>& p
    ) const
    {
      const auto other = p.on_frame(get_frame());
      return vector_t<decltype(value_type{} - U{})>{
        this->x() - other.x(),
        this->y() - other.y(),
        this->z() - other.z(),
        get_frame()
      };
    }

    // * point-vector sum/subtraction produces a point

    template <std::common_with<value_type> U>
    auto operator+(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return point_t<decltype(value_type{} + U{})>{
        this->x() + other.x(),
        this->y() + other.y(),
        this->z() + other.z(),
        get_frame()
      };
    }

    template <std::common_with<value_type> U>
    auto operator-(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return point_t<decltype(value_type{} - U{})>{
        this->x() - other.x(),
        this->y() - other.y(),
        this->z() - other.z(),
        get_frame()
      };
    }

    template <std::convertible_to<value_type> U>
    point_t& operator+=(
      const vector_t<U>& v
    )
    {
      const auto other = v.on_frame(get_frame());
      this->x() += other.x();
      this->y() += other.y();
      this->z() += other.z();
      return *this;
    }

    template <std::convertible_to<value_type> U>
    point_t& operator-=(
      const vector_t<U>& v
    )
    {
      const auto other = v.on_frame(get_frame());
      this->x() -= other.x();
      this->y() -= other.y();
      this->z() -= other.z();
      return *this;
    }

    // * point-point comparison
    
    template <std::common_with<value_type> U>
    auto operator==(
      const point_t<U>& p
    ) const
    {
      const auto other = p.on_frame(get_frame());
      return
        other.x() == this->x() &&
        other.y() == this->y() &&
        other.z() == this->z();
    }

  };

}

#endif