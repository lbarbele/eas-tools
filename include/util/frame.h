#ifndef _util_frame_h
#define _util_frame_h

#include <memory>
#include <type_traits>

#include <util/matrix.h>
#include <util/rotation_matrix.h>

namespace util {

  class frame;
  using frame_ptr = std::shared_ptr<frame>;

  class frame {
  private:
    square_matrix_d<3> m_from;
    square_matrix_d<3> m_to;

    constexpr frame() {}
    constexpr frame(const frame& other) {}
    constexpr frame(frame&& other) {}

  public:

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

    // * create frame from sequence of rotations
    template <class... Args,
      typename = std::enable_if_t<(sizeof...(Args) >= 2) && (sizeof...(Args)%2 == 0)>
    >
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

  };

}

#endif