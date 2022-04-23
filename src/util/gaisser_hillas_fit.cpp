#include <util/gaisser_hillas_fit.h>

#include <cmath>

namespace util {

  gaisser_hillas_fit::gaisser_hillas_fit(
    const double nmax,
    const double x0,
    const double xmax,
    const double p1,
    const double p2,
    const double p3
  ) :
    m_param{nmax, x0, xmax, p1, p2, p3}
  {
  }

  double
  gaisser_hillas_fit::eval(
    const double x
  ) const
  {
    auto& [nmax, x0, xmax, p1, p2, p3] = m_param;

    const double y = x - x0;

    if (y <= 0) {
      return 0;
    }

    const double ymax = xmax - x0;
    const double dy = x - xmax;
    const double lambda = p1 + x*(p2 + x*p3);

    return nmax * std::pow(y/ymax, ymax/lambda) * std::exp(-dy/lambda);
  }

} // namespace util