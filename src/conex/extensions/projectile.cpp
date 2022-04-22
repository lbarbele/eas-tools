#include "projectile.h"

#include <cmath>

namespace conex::extensions{

  util::vector
  projectile::get_position()
  const
  {
    const double d = height + 6371315.0;
    const double z = std::sqrt(d*d - x*x - y*y) - 6371315.0;
    return {x, y, z};
  }

  util::vector
  projectile::get_momentum()
  const
  {
    double r = std::sqrt(x*x + y*y);
    if (r > 1e-20) { 
      double sinp = y/r;
      double cosp = x/r;
      double sint = r/(height + 6371315.0);
      double cost = std::sqrt((1+sint)*(1-sint));
      return util::vector::to_obs({Px, Py, Pz}, sinp, cosp, sint, cost);
    } else {
      return util::vector::to_obs({Px, Py, Pz}, 0, 1, 0, 1);
    }
  }

} // namespace conex::extensions