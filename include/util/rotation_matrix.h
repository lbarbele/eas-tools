#ifndef _util_rotation_matrix_h
#define _util_rotation_matrix_h

#include <cmath>

#include <util/matrix.h>
#include <util/units.h>
#include <util/math.h>

namespace util {

  // * definition of the rotation axes
  enum class axis : size_t {
    x = 0,
    y = 1,
    z = 2
  };

  // * rotation matrix class
  class rotation_matrix : public square_matrix_d<3> {
  public:

    // * build rotation matrix from single rotation around specified axis
    rotation_matrix(
      const axis ax,
      const units::radian_t<double> angle
    ) :
      square_matrix_d<3>(0)
    {
      // * sine/cosine of the rotation angle
      const double sin_angle = math::sin(angle);
      const double cos_angle = math::cos(angle);

      // * some helper indices
      const size_t iaxis = static_cast<size_t>(ax);
      const size_t i = (iaxis+1)%3;
      const size_t j = (iaxis+2)%3;
    
      // * diagonal elements
      at(iaxis, iaxis) = 1;
      at(i, i) = cos_angle;
      at(j, j) = cos_angle;

      // * off-diagonal elements
      at(i, j) = sin_angle;
      at(j, i) = -sin_angle;
    }
  };

  // * shorthand notation to create a rotation matrix
  rotation_matrix
  operator,(
    const axis ax,
    const concepts::quantity_compatible<units::radian_t<double>> auto& angle
  )
  {
    using namespace units::literals;
    return rotation_matrix(ax, angle);
  }

}

#endif