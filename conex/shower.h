#ifndef _conex_shower_h
#define _conex_shower_h

#include <cmath>
#include <vector>

#include "../util/constants.h"
#include "../util/gaisser_hillas_fit.h"

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
    util::gaisser_hillas_fit get_fit() const
    {
      util::gaisser_hillas_fit fit(Nmax, X0, Xmax, p1, p2, p3);
      fit.set_ndof(nX);
      fit.set_chi2(chi2);
      return fit;
    }

    // primary energy
    float get_lge() const
    {return lgE;}
    
    double get_energy_gev() const
    {return std::pow(10, lgE - 9);}

    // zenith/azimuth of shower axis
    float get_zenith_deg() const
    {return zenith;}
    double get_zenith_rad() const
    {return zenith*util::constants::pi/180.;}

    float get_azimuth_deg() const
    {return azimuth;}
    double get_azimuth_rad() const
    {return azimuth*util::constants::pi/180.;}

    // random seeds
    int get_seed2() const
    {return Seed2;}

    int get_seed3() const
    {return Seed3;}

    // data from the first interaction
    float get_first_interaction_depth() const
    {return Xfirst;}

    float get_first_interaction_height() const
    {return Hfirst;}

    float get_first_interaction_inelasticty() const
    {return XfirstIn;}

    // impact parameter
    double get_altitude() const
    {return altitude;}

    // gaisser hillas fit (6 parameters)
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

    // real xmax
    float get_xmx() const
    {return Xmx;}

    // real nmax
    float get_nmx() const
    {return Nmx;}

    // real xmax of dEdX profile
    float get_xmx_dedx() const
    {return XmxdEdX;}

    // real max value of dEdX profile
    float get_dedx_mx() const
    {return dEdXmx;}

    // cpu time in seconds
    float get_cpu_time() const
    {return cpuTime;}

    // number of points in the longitudinal profiles
    int get_nx() const
    {return nX;}

    // actual profiles
    const float* get_depths() const
    {return X;}

    const float* get_charged() const
    {return N;}

    const float* get_heights() const
    {return H;}

    const float* ged_distances() const
    {return D;}

    const float* get_dedx() const
    {return dEdX;}

    const float* get_muons() const
    {return Mu;}

    const float* get_photons() const
    {return Gamma;}

    const float* get_electrons() const
    {return Electrons;}

    const float* get_hadrons() const
    {return Hadrons;}

    const float* get_mpd() const
    {return dMu;}

    // energy at ground
    float get_ground_energy_em() const
    {return EGround[0];}

    float get_ground_energy_hadrons() const
    {return EGround[1];}

    float get_ground_energy_muons() const
    {return EGround[2];}
  };

} // namespace conex

#endif // _conex_shower_h