#include <cmath>

#include <util/math.h>

namespace util::math {

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

  double
  gaisser_hillas(
    const double* x,
    const double* p
  )
  {
    return gaisser_hillas(*x, p[0], p[1], p[2], p[3]);
  }

  double
  double_gaisser_hillas(
    const double x,
    const double ecal,
    const double w,
    const double x01,
    const double xmax1,
    const double lambda1,
    const double x02,
    const double xmax2,
    const double lambda2
  )
  {
    return w*gaisser_hillas(x, ecal, x01, xmax1, lambda1) +
      (1-w)*gaisser_hillas(x, ecal, x02, xmax2, lambda2);
  }

  double
  double_gaisser_hillas(
    const double* x,
    const double* p
  )
  {
    return double_gaisser_hillas(*x, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
  }

  double
  gaisser_hillas_sixpar(
    const double x,
    const double ymax,
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
    return ymax * std::pow((x-x0)/(xmax-x0), (xmax-x0)/l) * std::exp((xmax-x)/l);
  }

  double
  gaisser_hillas_sixpar(
    const double* x,
    const double* p
  )
  {
    return gaisser_hillas_sixpar(*x, p[0], p[1], p[2], p[3], p[4], p[5]);
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