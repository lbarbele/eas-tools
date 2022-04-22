#ifndef _util_gaisser_hillas_fit_h
#define _util_gaisser_hillas_fit_h

namespace util {

  class gaisser_hillas_fit {
  private:
    double m_nmax;
    double m_x0;
    double m_xmax;
    double m_p1;
    double m_p2;
    double m_p3;
    double m_dev;
    double m_chi2;
    int m_ndof;

    void update_lambda(const double p1, const double p2, const double p3);

  public:
    gaisser_hillas_fit(const double nmax, const double x0, const double xmax, const double p1, const double p2 = 0, const double p3 = 0);

    // evaluators
    double eval(const double x) const;

    double operator()(const double x) const
    {return eval(x);}

    // getters/setters
    void set_nmax(const double nmax)
    {m_nmax = nmax;}

    void set_x0(const double x0)
    {m_x0 = x0;}

    void set_xmax(const double xmax)
    {m_xmax = xmax;}

    void set_lambda(const double lambda)
    {m_p1 = lambda;}

    void set_p1(const double p1)
    {m_p1 = p1;}

    void set_p2(const double p2)
    {m_p2 = p2;}

    void set_p3(const double p3)
    {m_p3 = p3;}

    void set_dev(const double dev)
    {m_dev = dev;}

    void set_chi2(const double chi2)
    {m_chi2 = chi2;}

    void set_ndof(const int ndof)
    {m_ndof = ndof;}

    // getters
    double get_nmax() const
    {return m_nmax;}

    double get_x0() const
    {return m_x0;}

    double get_xmax() const
    {return m_xmax;}

    double get_p1() const
    {return m_p1;}

    double get_p2() const
    {return m_p2;}

    double get_p3() const
    {return m_p3;}

    double get_dev() const
    {return m_dev;}

    double get_chi2() const
    {return m_chi2;}

    int get_ndof() const
    {return m_ndof;}

  };

} // namespace util

#endif