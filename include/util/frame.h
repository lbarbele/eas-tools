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
    square_matrix_d<3> m_from{1, 0, 0, 0, 1, 0, 0, 0, 1};
    square_matrix_d<3> m_to{1, 0, 0, 0, 1, 0, 0, 0, 1};

    // * coordinates of origin wrt standard frame
    coordinates_t<double> m_origin{0, 0, 0};

    constexpr frame() {}
    constexpr frame(const frame& other) = delete;
    constexpr frame(frame&& other) = delete;

    void set_rotation_matrix(
      const square_matrix_d<3>& rot_to,
      const frame_ptr& base_frame
    )
    {
      m_to = rot_to*base_frame->to();
      m_from = m_to.transpose();
    }

    void set_origin_coordinates(
      const coordinates_t<double>& new_origin,
      const frame_ptr& base_frame
    )
    {
      m_origin = base_frame->from() * new_origin;
      m_origin.x() += base_frame->origin().x();
      m_origin.y() += base_frame->origin().y();
      m_origin.z() += base_frame->origin().z();
    }

  public:

    // - Static creation methods

    // * create the default frame
    static frame_ptr create()
    {
      return frame_ptr(new frame);
    }

    // * create copy of frame
    static frame_ptr create(
      const frame_ptr& base_frame
    )
    {
      return base_frame;
    }

    // * create displaced frame, with same orientation
    static frame_ptr create(
      const coordinates_t<double>& new_origin,
      const frame_ptr& base_frame
    )
    {
      frame_ptr f = create();
      f->m_from = base_frame->from();
      f->m_to = base_frame->to();
      f->set_origin_coordinates(new_origin, base_frame);
      return f;
    }

    // * create rotated frame, with same origin
    static frame_ptr create(
      const square_matrix_d<3>& rot_to,
      const frame_ptr& base_frame
    )
    {
      frame_ptr f = create();
      f->set_rotation_matrix(rot_to, base_frame);
      f->m_origin = base_frame->origin();
      return f;
    }

    // * create frame from rotation and displacement
    static frame_ptr create(
      const square_matrix_d<3>& rot_to,
      const coordinates_t<double>& new_origin,
      const frame_ptr& base_frame
    )
    {
      frame_ptr f = create();
      f->set_rotation_matrix(rot_to, base_frame);
      f->set_origin_coordinates(new_origin, base_frame);
      return f;
    }

    static frame_ptr create(
      const coordinates_t<double>& new_origin,
      const square_matrix_d<3>& rot_to,
      const frame_ptr& base_frame
    )
    {
      return create(rot_to, new_origin, base_frame);
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