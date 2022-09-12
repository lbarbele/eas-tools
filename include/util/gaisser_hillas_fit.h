#ifndef _util_gaisser_hillas_fit_h
#define _util_gaisser_hillas_fit_h

#include <array>
#include <algorithm>
#include <cmath>

namespace util {

  class gaisser_hillas_fit {
  private:
    std::array<double, 6> m_param;

  public:
    gaisser_hillas_fit(
      const double nmax = 0,
      const double x0 = 0,
      const double xmax = 0,
      const double p1 = 0,
      const double p2 = 0,
      const double p3 = 0
    ) :
      m_param{nmax, x0, xmax, p1, p2, p3}
    {}

    template<class array_t>
    gaisser_hillas_fit(const array_t& v)
    {set_parameters(v);}

    // evaluators
    double
    eval(
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

    double operator()(const double x) const
    {return eval(x);}

    // set all parameters 
    template<class array_t>
    void set_parameters(const array_t& v)
    {
      std::copy(v.begin(), v.end(), m_param.begin());
    }

    // getters/setters for single parameters
    void set_nmax(const double nmax)
    {m_param[0] = nmax;}
    double get_nmax() const
    {return m_param[0];}

    void set_x0(const double x0)
    {m_param[1] = x0;}
    double get_x0() const
    {return m_param[1];}

    void set_xmax(const double xmax)
    {m_param[2] = xmax;}
    double get_xmax() const
    {return m_param[2];}

    void set_lambda(const double lambda)
    {m_param[3] = lambda;}
    double get_lambda() const
    {return m_param[3];}

    void set_p1(const double p1)
    {m_param[4] = p1;}
    double get_p1() const
    {return m_param[4];}

    void set_p2(const double p2)
    {m_param[5] = p2;}
    double get_p2() const
    {return m_param[5];}

    void set_p3(const double p3)
    {m_param[6] = p3;}
    double get_p3() const
    {return m_param[6];}

  };

} // namespace util

#endif