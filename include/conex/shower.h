#ifndef _conex_shower_h
#define _conex_shower_h

#include <TGraph.h>

#include <util/constants.h>
#include <util/frame.h>
#include <util/gaisser_hillas_fit.h>
#include <util/math.h>
#include <util/units.h>
#include <util/vector.h>

namespace conex {

  struct shower {
    friend class file;

  private:
    static const unsigned int max_entries = 5000;
    
    float                                         lgE;
    units::degree_t<float>                        zenith;
    units::degree_t<float>                        azimuth;
    int                                           Seed2;
    int                                           Seed3;
    units::grams_per_squared_centimeter_t<float>  Xfirst;
    units::meter_t<float>                         Hfirst;
    float                                         XfirstIn;
    units::meter_t<double>                        altitude;
    units::grams_per_squared_centimeter_t<float>  X0;
    units::grams_per_squared_centimeter_t<float>  Xmax;
    float                                         Nmax;
    units::grams_per_squared_centimeter_t<float>  p1;
    float                                         p2;
    units::squared_centimeters_per_gram_t<float>  p3;
    float                                         chi2;
    units::grams_per_squared_centimeter_t<float>  Xmx;
    float                                         Nmx;
    units::grams_per_squared_centimeter_t<float>  XmxdEdX;
    units::gev_per_gcm_t<float>                   dEdXmx;
    units::second_t<float>                        cpuTime;
    int                                           nX;

    units::grams_per_squared_centimeter_t<float>  X[max_entries];
    units::meter_t<float>  H[max_entries];
    units::meter_t<float>  D[max_entries];
    units::gev_per_gcm_t<float>  dEdX[max_entries];
    float  N[max_entries];
    float  Mu[max_entries];
    float  Gamma[max_entries];
    float  Electrons[max_entries];
    float  Hadrons[max_entries];
    float  dMu[max_entries];

    units::gigaelectron_volt_t<float>  EGround[3];

  public:

    // * direct access to CONEX scalar data

    // primary energy
    units::energy_t get_energy() const
    {return units::gigaelectron_volt_t<double>(util::math::pow(10, lgE - 9));}

    // zenith/azimuth of shower axis
    units::angle_t get_zenith() const
    {return zenith;}

    units::angle_t get_azimuth() const
    {return azimuth;}

    // random seeds
    int get_seed2() const
    {return Seed2;}

    int get_seed3() const
    {return Seed3;}

    // data from the first interaction
    units::depth_t get_first_interaction_depth() const
    {return Xfirst;}

    units::length_t get_first_interaction_height() const
    {return Hfirst;}

    double get_first_interaction_inelasticty() const
    {return XfirstIn;}

    // impact parameter
    units::length_t get_altitude() const
    {return altitude;}

    // gaisser hillas fit (6 parameters)
    units::depth_t get_x0() const
    {return X0;}

    units::depth_t get_xmax() const
    {return Xmax;}

    double get_nmax() const
    {return Nmax;}

    units::depth_t get_p1() const
    {return p1;}

    double get_p2() const
    {return p2;}

    units::squared_centimeters_per_gram_t<float> get_p3() const
    {return p3;}

    double get_chi2() const
    {return chi2;}

    // real xmax
    units::depth_t get_xmx() const
    {return Xmx;}

    // real nmax
    double get_nmx() const
    {return Nmx;}

    // real xmax of dEdX profile
    units::depth_t get_xmx_dedx() const
    {return XmxdEdX;}

    // real max value of dEdX profile
    units::energy_deposit_t get_dedx_mx() const
    {return dEdXmx;}

    // cpu time in seconds
    units::time_t get_cpu_time() const
    {return cpuTime;}

    // number of points in the longitudinal profiles
    int get_nx() const
    {return nX;}

    // * profiles

    // (slant) atmospheric depth profile
    const units::grams_per_squared_centimeter_t<float>* get_depths() const
    {return X;}

    // height (above sea level) profile
    const units::meter_t<float>* get_heights() const
    {return H;}

    // distances to impact point
    const units::meter_t<float>* get_distances() const
    {return D;}

    // energy deposit profile
    const units::gev_per_gcm_t<float>* get_dedx() const
    {return dEdX;}

    // number of charged particles
    const float* get_charged() const
    {return N;}

    // number of muons
    const float* get_muons() const
    {return Mu;}

    // number of photons
    const float* get_photons() const
    {return Gamma;}

    // number of electrons + positrons
    const float* get_electrons() const
    {return Electrons;}

    // number of hadrons
    const float* get_hadrons() const
    {return Hadrons;}

    // muon production depth profile
    const float* get_mpd() const
    {return dMu;}

    // * energy at ground
    
    units::energy_t get_ground_energy_em() const
    {return EGround[0];}

    units::energy_t get_ground_energy_hadrons() const
    {return EGround[1];}

    units::energy_t get_ground_energy_muons() const
    {return EGround[2];}

    // * access to the gaisser hillas fit

    util::gaisser_hillas_fit get_fit() const
    {
      using namespace units::literals;
      return util::gaisser_hillas_fit(Nmax, X0/1_gcm2, Xmax/1_gcm2, p1/1_gcm2, p2, p3*1_gcm2);
    }

    // * generate profile TGraphs

    TGraph
    graph_dedx()
    const
    {
      using namespace units::literals;
      TGraph g(get_nx() - 1);
      for (int i = 0; i < get_nx() - 1; ++i) {
        const units::depth_t x = 0.5*(get_depths()[i] + get_depths()[i+1]);
        g.SetPoint(i, x/1_gcm2, get_dedx()[i] * (1_gcm2/1_GeV));
      }
      return g;
    }

    // * get vector holding the shower axis

    util::vector_d
    get_axis()
    const
    {
      const auto theta = get_zenith();
      const auto phi = get_azimuth();
      return util::vector_d(
        util::math::sin(theta) * util::math::cos(phi),
        util::math::sin(theta) * util::math::sin(phi),
        util::math::cos(theta),
        util::frame::conex_observer
      );
    }
  };

} // namespace conex

#endif // _conex_shower_h