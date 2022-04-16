#ifndef _conex_shower_h
#define _conex_shower_h

namespace conex {

  struct shower {
    friend class file;

  private:
    static const unsigned int max_entries = 5000;
    
    float  lgE;
    float  zenith;
    float  azimuth;
    int    Seed2;
    int    Seed3;
    float  Xfirst;
    float  Hfirst;
    float  XfirstIn;
    double altitude;
    float  X0;
    float  Xmax;
    float  Nmax;
    float  p1;
    float  p2;
    float  p3;
    float  chi2;
    float  Xmx;
    float  Nmx;
    float  XmxdEdX;
    float  dEdXmx;
    float  cpuTime;
    int    nX;

    float  X[max_entries];
    float  N[max_entries];
    float  H[max_entries];
    float  D[max_entries];
    float  dEdX[max_entries];
    float  Mu[max_entries];
    float  Gamma[max_entries];
    float  Electrons[max_entries];
    float  Hadrons[max_entries];
    float  dMu[max_entries];
    float  EGround[3];

  public:
    float get_lge() const
    {return lgE;}

    float get_zenith() const
    {return zenith;}

    float get_azimuth() const
    {return azimuth;}

    int get_seed2() const
    {return Seed2;}

    int get_seed3() const
    {return Seed3;}

    float get_first_x() const
    {return Xfirst;}

    float get_first_h() const
    {return Hfirst;}

    float get_first_inel() const
    {return XfirstIn;}

    double get_altitude() const
    {return altitude;}

    float get_x0() const
    {return X0;}

    float get_xmax() const
    {return Xmax;}

    float get_nmax() const
    {return Nmax;}

    float get_p1() const
    {return p1;}

    float get_p2() const
    {return p2;}

    float get_p3() const
    {return p3;}

    float get_chi2() const
    {return chi2;}

    float get_xmx() const
    {return Xmx;}

    float get_nmx() const
    {return Nmx;}

    float get_xmx_dedx() const
    {return XmxdEdX;}

    float get_dedx_mx() const
    {return dEdXmx;}

    float get_cpu_time() const
    {return cpuTime;}

    int get_nx() const
    {return nX;}


  };

} // namespace conex

#endif // _conex_shower_h