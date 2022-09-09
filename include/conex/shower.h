#ifndef _conex_shower_h
#define _conex_shower_h

#include <cmath>
#include <vector>

#include <TGraph.h>

#include <util/constants.h>
#include <util/gaisser_hillas_fit.h>
#include <util/units.h>
#include <util/vector.h>

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

    TGraph graph_dedx() const;

    util::vector_d get_axis() const;

    util::gaisser_hillas_fit get_fit() const
    {return util::gaisser_hillas_fit(Nmax, X0, Xmax, p1, p2, p3);}

    // primary energy
    auto get_energy() const
    {return units::gigaelectron_volt_t<double>(std::pow(10, lgE - 9));}

    // zenith/azimuth of shower axis
    auto get_zenith() const
    {return units::degree_t<double>(zenith);}

    auto get_azimuth() const
    {return units::degree_t<double>(azimuth);}

    // random seeds
    int get_seed2() const
    {return Seed2;}

    int get_seed3() const
    {return Seed3;}

    // data from the first interaction
    auto get_first_interaction_depth() const
    {return units::grams_per_squared_centimeter_t<double>(Xfirst);}

    auto get_first_interaction_height() const
    {return units::meter_t<double>(Hfirst);}

    float get_first_interaction_inelasticty() const
    {return XfirstIn;}

    // impact parameter
    auto get_altitude() const
    {return units::meter_t<double>(altitude);}

    // gaisser hillas fit (6 parameters)
    auto get_x0() const
    {return units::grams_per_squared_centimeter_t<double>(X0);}

    auto get_xmax() const
    {return units::grams_per_squared_centimeter_t<double>(Xmax);}

    double get_nmax() const
    {return Nmax;}

    auto get_p1() const
    {return units::grams_per_squared_centimeter_t<double>(p1);}

    double get_p2() const
    {return p2;}

    auto get_p3() const
    {return p3 / units::grams_per_squared_centimeter_t<double>(1);}

    float get_chi2() const
    {return chi2;}

    // real xmax
    auto get_xmx() const
    {return units::grams_per_squared_centimeter_t<double>(Xmx);}

    // real nmax
    float get_nmx() const
    {return Nmx;}

    // real xmax of dEdX profile
    auto get_xmx_dedx() const
    {return units::grams_per_squared_centimeter_t<double>(XmxdEdX);}

    // real max value of dEdX profile
    auto get_dedx_mx() const
    {
      using namespace units::literals;
      return dEdXmx * (1_GeV/1_gcm2);
    }

    // cpu time in seconds
    auto get_cpu_time() const
    {return units::second_t<double>(cpuTime);}

    // number of points in the longitudinal profiles
    int get_nx() const
    {return nX;}

    // actual profiles
    auto get_depths() const
    {
      using quantity_type = units::grams_per_squared_centimeter_t<float>;
      return (quantity_type*)X;
    }

    const float* get_charged() const
    {return N;}

    auto get_heights() const
    {return (units::meter_t<float>*)H;}

    auto get_distances() const
    {return (units::meter_t<float>*)D;}

    auto get_dedx() const
    {
      using namespace units;
      using unit_type = make_unit<gigaelectron_volt, inverse<grams_per_squared_centimeter>>;
      using quantity_type = quantity<unit_type, float>;
      return (quantity_type*)dEdX;
    }

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
    auto get_ground_energy_em() const
    {return units::gigaelectron_volt_t<double>(EGround[0]);}

    auto get_ground_energy_hadrons() const
    {return units::gigaelectron_volt_t<double>(EGround[1]);}

    auto get_ground_energy_muons() const
    {return units::gigaelectron_volt_t<double>(EGround[2]);}
  };

} // namespace conex

#endif // _conex_shower_h