#ifndef _util_frame_h
#define _util_frame_h

#include <memory>
#include <type_traits>

#include <util/matrix.h>
#include <util/rotation_matrix.h>
#include <util/constants.h>

namespace util {

  class frame;
  using frame_ptr = std::shared_ptr<frame>;

  class frame {
  private:
    square_matrix_d<3> m_from;
    square_matrix_d<3> m_to;

    constexpr frame() : m_from(0), m_to(0) {}
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

    // * create frame from a single rotation around specified axis
    static frame_ptr create(
      const axis ax,
      const double angle
    )
    {
      frame_ptr f(new frame());
      f->m_to = rotation_matrix(ax, angle);
      f->m_from = f->m_to.transpose();
      return f;
    }

    // * create frame from a sequence of rotations
    template <class... Args, std::enable_if_t<(sizeof...(Args) >= 2) && (sizeof...(Args)%2 == 0), bool> = true >
    static frame_ptr create(
      const axis ax,
      const double angle,
      Args... args
    )
    {
      auto f = create(args...);
      f->m_to *= rotation_matrix(ax, angle);
      return f;
    }

    // - Methods

    // * access rotation matrices
    constexpr const square_matrix_d<3>& to() const
    {return m_to;}

    constexpr const square_matrix_d<3>& from() const
    {return m_from;}

    // - Default frames
    const inline static frame_ptr corsika_observer = create();
    const inline static frame_ptr conex_observer = create(axis::z, -constants::pi/2.0);
    const inline static frame_ptr standard = corsika_observer;

  };

}

#endif