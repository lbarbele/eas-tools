#include <cmath>

#include <util/math.h>

namespace util::math {

  double
  gaisser_hillas(
    const double x,
    const double nmax,
    const double x0,
    const double xmax,
    const double p1,
    const double p2,
    const double p3
  )
  {
    if (x <= x0) {
      return 0;
    }

    const double l = p1 + x*(p2 + x*p3);
    const double a = (xmax - x0) / l;
    const double r = (x - x0) / (xmax - x0);
    return nmax * std::exp(a * (1 - r + std::log(r)));
  }

  double
  gaisser_hillas_ecal(
    const double x,
    const double ecal,
    const double x0,
    const double xmax,
    const double l
  )
  {
    if (x <= x0) {
      return 0;
    }
    
    const double z = (x - x0) / l;
    const double a = (xmax - x0) / l;
    return (ecal/l) * std::exp(-std::lgamma(a+1) - z + a*std::log(z));
  }

  double
  usp_ecal(
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

} // namespace util::math