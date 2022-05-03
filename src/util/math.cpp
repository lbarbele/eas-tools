#include <cmath>

#include <util/math.h>

namespace util::math {

  // gaisser-hillas function: explicit parameters
  double
  gaisser_hillas(
    const double x,
    const double ecal,
    const double x0,
    const double xmax,
    const double l
  )
  {
    const double z = (x - x0)/l;
    const double am1 = (xmax - x0)/l;
    return z <= 0? 0 : ecal * std::pow(z, am1) * std::exp(-z) / (l * std::tgamma(1+am1));
  }

  // gaisser-hillas function: compatible with ROOT's TF1
  double
  gaisser_hillas(
    const double* x,
    const double* p
  )
  {
    return gaisser_hillas(*x, p[0], p[1], p[2], p[3]);
  }

  // universal shower profile function: explicit parameters
  double
  usp_function(
    const double x,
    const double ecal,
    const double xmax,
    const double l,
    const double r
  )
  {
    const double invr2 = 1.0/(r*r);
    const double z = (1.0/r + (x-xmax)/l) / r;
    return z<=0? 0 : ecal*std::pow(z,invr2)*std::exp(-z)/(l*r*std::tgamma(1+invr2));
  }

  // universal shower profile function: compatible with ROOT's TF1
  double
  usp_function(
    const double* x,
    const double* p
  )
  {
    return usp_function(*x, p[0], p[1], p[2], p[3]);
  }

} // namespace util::math