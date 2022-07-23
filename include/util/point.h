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
        // rotate to the standard frame
        *this = get_frame()->from() * (*this);
        // move origin to the new frame
        this->x() += get_frame()->origin()[0] - frame->origin()[0];
        this->y() += get_frame()->origin()[1] - frame->origin()[1];
        this->z() += get_frame()->origin()[2] - frame->origin()[2];
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

    // * distance to another point
    template <class U>
    decltype(T{} - U{}) distance(point_t<U> p) const
    {
      p.on_frame(get_frame());
      return std::hypot(this->x()-p.x(), this->y()-p.y(), this->z()-p.z());
    }

    // * point-point subtraction produces a vector
    template <class U>
    vector_t<decltype(T{}-U{})>
    operator-(
      point_t<U> p
    ) const
    {
      p.set_frame(get_frame());
      return {this->x() - p.x(), this->y() - p.y(), this->z() - p.z(), get_frame()};
    }

    // * point-vector sum/subtraction produces a point
    template <class U>
    point_t<decltype(T{}+U{})>
    operator+(
      vector_t<U> v
    ) const
    {
      v.set_frame(get_frame());
      return {this->x() + v[0], this->y() + v[1], this->z() + v[2], get_frame()};
    }

    template <class U>
    point_t<T>&
    operator+=(
      vector_t<U> v
    )
    {
      v.set_frame(get_frame());
      this->x() += v[0];
      this->y() += v[1];
      this->z() += v[2];
      return *this;
    }

    template <class U>
    point_t<decltype(T{}-U{})>
    operator-(
      vector_t<U> v
    ) const
    {
      v.set_frame(get_frame());
      return {this->x() - v[0], this->y() - v[1], this->z() - v[2], get_frame()};
    }

    template <class U>
    point_t<T>&
    operator-=(
      vector_t<U> v
    )
    {
      v.set_frame(get_frame());
      this->x() -= v[0];
      this->y() -= v[1];
      this->z() -= v[2];
      return *this;
    }

  };

  // - Aliases
  using point_d = point_t<double>;
  using point_f = point_t<float>;

}

#endif