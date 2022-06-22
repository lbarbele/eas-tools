#ifndef _util_rotation_matrix_h
#define _util_rotation_matrix_h

#include <cmath>

#include <util/matrix.h>

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
    rotation_matrix(
      const axis ax,
      const double angle
    ) :
      square_matrix_d<3>(0)
    {
      // * sine/cosine of the rotation angle
      const double sin_angle = std::sin(angle);
      const double cos_angle = std::cos(angle);

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

}

#endif