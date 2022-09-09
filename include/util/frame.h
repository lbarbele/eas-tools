#ifndef _util_frame_h
#define _util_frame_h

#include <array>
#include <concepts>
#include <memory>
#include <type_traits>

#include <any>

#include <util/constants.h>
#include <util/coordinates.h>
#include <util/matrix.h>
#include <util/rotation_matrix.h>

namespace util {

  // - forward declaration and aliases

  // * actual frame definition is hidden in _impl namespace
  namespace _impl {
    template <class T> class frame;
  }

  // * frame_scale is a wrapper for the frame scale type providing the necessary type decays
  template <class T>
  struct frame_scale;

  // * alias for the decayed type of the frame scale
  template <class T>
  using frame_scale_t = typename frame_scale<T>::type;

  // * import the frame to the util NS, but with the scale type decayed
  template <class T>
  using frame = _impl::frame<frame_scale_t<T>>;

  // * smart pointer for a frame
  template <class T>
  using frame_ptr = typename frame<T>::ptr_type;

  // - implementation scale type wrapper

  // * specialization for arithmetic types
  template <units::concepts::arithmetic T>
  struct frame_scale<T> {
    using type = std::remove_cvref_t<T>;
  };

  // * specialization for quantities with dimensions
  template <units::concepts::quantity Q>
  struct frame_scale<Q> {
    using value_type = typename Q::value_type;
    using unit_type = units::make_unit<typename Q::unit_type, units::ratio_power<typename Q::unit_type::factor, -1>>;
    using type = units::quantity<unit_type, value_type>;
  };

  // - implementation of the frame class

  namespace _impl {

    template <class T>
    class frame {
    public:

      using ptr_type = std::shared_ptr<frame>;
      using scale = T;

    private:

      // * rotation matrices with respect to the standard frame
      square_matrix_d<3> m_from{1, 0, 0, 0, 1, 0, 0, 0, 1};
      square_matrix_d<3> m_to{1, 0, 0, 0, 1, 0, 0, 0, 1};

      // * coordinates of origin wrt standard frame
      coordinates_t<scale> m_origin{scale(0), scale(0), scale(0)};

      constexpr frame() {}
      constexpr frame(const frame& other) = delete;
      constexpr frame(frame&& other) = delete;

      void set_rotation_matrix(
        const square_matrix_d<3>& rot_to,
        const ptr_type& base_frame
      )
      {
        m_to = rot_to*base_frame->to();
        m_from = m_to.transpose();
      }

      void set_origin_coordinates(
        const coordinates_t<scale>& new_origin,
        const ptr_type& base_frame
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
      static ptr_type create()
      {return ptr_type(new frame);}

      // * create copy of frame
      static ptr_type create(const ptr_type& base_frame)
      {return base_frame;}

      // * create displaced frame, with same orientation
      static ptr_type create(
        const coordinates_t<scale>& new_origin,
        const ptr_type& base_frame
      )
      {
        ptr_type f = create();
        f->m_from = base_frame->from();
        f->m_to = base_frame->to();
        f->set_origin_coordinates(new_origin, base_frame);
        return f;
      }

      // * create rotated frame, with same origin
      static ptr_type create(
        const square_matrix_d<3>& rot_to,
        const ptr_type& base_frame
      )
      {
        ptr_type f = create();
        f->set_rotation_matrix(rot_to, base_frame);
        f->m_origin = base_frame->origin();
        return f;
      }

      // * create frame from rotation and displacement
      static ptr_type create(
        const square_matrix_d<3>& rot_to,
        const coordinates_t<scale>& new_origin,
        const ptr_type& base_frame
      )
      {
        ptr_type f = create();
        f->set_rotation_matrix(rot_to, base_frame);
        f->set_origin_coordinates(new_origin, base_frame);
        return f;
      }

      static ptr_type create(
        const coordinates_t<scale>& new_origin,
        const square_matrix_d<3>& rot_to,
        const ptr_type& base_frame
      )
      {return create(rot_to, new_origin, base_frame);}

      // - Methods

      // * access rotation matrices
      constexpr const square_matrix_d<3>& to() const
      {return m_to;}

      constexpr const square_matrix_d<3>& from() const
      {return m_from;}

      // * access coordinates of origin
      constexpr const coordinates_t<scale>& origin() const
      {return m_origin;}

      // - Default frames
      
      const inline static ptr_type standard = create();
      const inline static ptr_type corsika_observer = standard;
      const inline static ptr_type conex_observer = create((axis::z, -constants::pi/2.0), standard);

    };

  }

}

#endif