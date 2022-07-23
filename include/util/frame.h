#ifndef _util_frame_h
#define _util_frame_h

#include <array>
#include <memory>
#include <type_traits>

#include <util/constants.h>
#include <util/coordinates.h>
#include <util/matrix.h>
#include <util/rotation_matrix.h>

namespace util {

  class frame;
  using frame_ptr = std::shared_ptr<frame>;

  class frame {
  private:
    // * rotation matrices with respect to the standard frame
    square_matrix_d<3> m_from{0};
    square_matrix_d<3> m_to{0};

    // * coordinates of origin wrt standard frame
    coordinates_t<double> m_origin{0, 0, 0};

    constexpr frame() {}
    constexpr frame(const frame& other) {}
    constexpr frame(frame&& other) {}

  public:

    // * create the default frame
    static frame_ptr create()
    {
      frame_ptr f(new frame());
      f->m_from(0, 0) = f->m_from(1, 1) = f->m_from(2, 2) = 1;
      f->m_to(0, 0) = f->m_to(1, 1) = f->m_to(2, 2) = 1;
      return f;
    }

    // * create frame from a transformation matrix ("to" matrix) and a starting frame
    // * origin of frames coincide
    static frame_ptr create(
      const square_matrix_d<3>& rot_to,
      const frame_ptr& base_frame
    )
    {
      frame_ptr f(new frame());
      f->m_to = rot_to*base_frame->to();
      f->m_from = f->m_to.transpose();
      f->m_origin = base_frame->origin();
      return f;
    }

    // * create frame from displacement only
    // * orientations coincide
    static frame_ptr create(
      const coordinates_t<double>& new_origin,
      const frame_ptr& base_frame
    )
    {
      // create empty frame
      frame_ptr f = create();

      // copy rotation matrices from the base frame
      f->m_from = base_frame->from();
      f->m_to = base_frame->to();

      // set origin as base_frame origin + new_origin
      const auto new_origin_std = base_frame->from() * new_origin;
      f->m_origin.x() = base_frame->origin().x() + new_origin_std.x();
      f->m_origin.y() = base_frame->origin().y() + new_origin_std.y();
      f->m_origin.z() = base_frame->origin().z() + new_origin_std.z();

      return f;
    }

    // * create frame from rotation and displacement
    static frame_ptr create(
      const square_matrix_d<3>& rot_to,
      const coordinates_t<double>& new_origin,
      const frame_ptr& base_frame
    )
    {
      // create a rotated frame with same origin as the base frame
      frame_ptr f = create(rot_to, base_frame);

      // set origin as base_frame origin + new_origin
      const auto new_origin_std = base_frame->from() * new_origin;
      f->m_origin.x() += new_origin_std.x();
      f->m_origin.y() += new_origin_std.y();
      f->m_origin.z() += new_origin_std.z();

      return f;
    }

    // - Methods

    // * access rotation matrices
    constexpr const square_matrix_d<3>& to() const
    {return m_to;}

    constexpr const square_matrix_d<3>& from() const
    {return m_from;}

    // * access coordinates of origin
    constexpr const coordinates_t<double>& origin() const
    {return m_origin;}

    // - Default frames
    const inline static frame_ptr standard = create();
    const inline static frame_ptr corsika_observer = standard;
    const inline static frame_ptr conex_observer = create((axis::z, -constants::pi/2.0), standard);

  };

}

#endif