#include "gaisser_hillas_fit.h"

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
    m_nmax(nmax),
    m_x0(x0),
    m_xmax(xmax),
    m_p1(p1),
    m_p2(p2),
    m_p3(p3)
  {
  }

  double
  gaisser_hillas_fit::eval(
    const double x
  ) const
  {
    const double y = x - m_x0;

    if (y <= 0) {
      return 0;
    }

    const double ymax = m_xmax - m_x0;
    const double dy = x - m_xmax;
    const double lambda = m_p1 + x*(m_p2 + x*m_p3);

    return m_nmax * std::pow(y/ymax, ymax/lambda) * std::exp(-dy/lambda);
  }

} // namespace util