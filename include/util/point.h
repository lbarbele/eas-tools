#ifndef _util_point_h
#define _util_point_h

#include <array>
#include <cmath>

#include <util/coordinates.h>
#include <util/frame.h>
#include <util/vector.h>

namespace util {

  template <class T>
  class point_t : public point_base_t, public coordinates_t<T> {
  private:
    frame_ptr m_frame;

  public:

    // - Constructors

    // * undefined coordinates on given frame (defaults to standard frame)
    point_t(
      const frame_ptr frame = frame::standard
    ) :
      m_frame(frame)
    {}

    // * construct a point from a coordinates object and frame (explicit)
    point_t(
      const coordinates_t<T>& coordinates,
      const frame_ptr frame
    ) :
      coordinates_t<T>(coordinates),
      m_frame(frame)
    {}

    point_t(
      coordinates_t<T>&& coordinates,
      const frame_ptr frame
    ) :
      coordinates_t<T>(coordinates),
      m_frame(frame)
    {}

    // * construct a point with given coordinates and frame
    point_t(
      const T& x,
      const T& y,
      const T& z,
      const frame_ptr frame = frame::standard
    ) :
      point_t({x, y, z}, frame)
    {}

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
        // copy coordinates of this point but rotated to coincide with the std frame
        coordinates_t<T> new_coordinates = get_frame()->from() * (*this);
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
    point_t<T> on_frame(
      const frame_ptr& frame
    ) const
    {
      point_t other = (*this);
      other.set_frame(frame);
      return other;
    }

    // - Operations

    // * distance to another point
    template <class U, class R = decltype(T{} - U{})>
    R
    distance(
      const point_t<U>& p
    ) const
    {
      const auto other = p.on_frame(get_frame());
      return coordinates_t<R>(
        this->x() - other.x(),
        this->y() - other.y(),
        this->z() - other.z()
      ).get_r();
    }

    // * point-point subtraction produces a vector
    template <class U, class R = decltype(T{} - U{})>
    vector_t<R>
    operator-(
      const point_t<U>& p
    ) const
    {
      const auto other = p.on_frame(get_frame());
      return {
        this->x() - other.x(),
        this->y() - other.y(),
        this->z() - other.z(),
        get_frame()
      };
    }

    // * point-vector sum/subtraction produces a point
    template <class U, class R = decltype(T{} + U{})>
    point_t<R>
    operator+(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return {
        this->x() + other.x(),
        this->y() + other.y(),
        this->z() + other.z(),
        get_frame()
      };
    }

    template <class U, class R = decltype(T{} - U{})>
    point_t<R>
    operator-(
      const vector_t<U>& v
    ) const
    {
      const auto other = v.on_frame(get_frame());
      return {
        this->x() - other.x(),
        this->y() - other.y(),
        this->z() - other.z(),
        get_frame()
      };
    }

    template <class U>
    point_t<T>&
    operator+=(
      const vector_t<U>& v
    )
    {
      const auto other = v.on_frame(get_frame());
      this->x() += other.x();
      this->y() += other.y();
      this->z() += other.z();
      return *this;
    }

    template <class U>
    point_t<T>&
    operator-=(
      const vector_t<U>& v
    )
    {
      const auto other = v.on_frame(get_frame());
      this->x() -= other.x();
      this->y() -= other.y();
      this->z() -= other.z();
      return *this;
    }

  };

  // - Aliases
  using point_d = point_t<double>;
  using point_f = point_t<float>;

}

#endif