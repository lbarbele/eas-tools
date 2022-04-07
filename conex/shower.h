#ifndef _conex_shower_h
#define _conex_shower_h

namespace conex {

  struct shower {
  private:
    static const unsigned int max_entries = 5000;
  public:
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
  };

} // namespace conex

#endif // _conex_shower_h