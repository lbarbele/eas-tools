#ifndef _conex_header_h
#define _conex_header_h

#include <array>
#include <iostream>

#include <util/math.h>
#include <util/units.h>

namespace conex {

  class header {
    friend class file;

  private:
    int    Seed1;
    int    Particle;
    double Alpha;
    double lgEmin;
    double lgEmax;
    double zMin;
    double zMax;
    float  Version;
    float  OutputVersion;
    int    HEModel;
    int    LEModel;
    float  HiLowEgy;
    float  hadCut;
    float  emCut;
    float  hadThr;
    float  muThr;
    float  emThr;
    float  haCut;
    float  muCut;
    float  elCut;
    float  gaCut;
    
    std::array<double, 31> lambdaLgE;
    std::array<double, 31> lambdaProton;
    std::array<double, 31> lambdaPion;
    std::array<double, 31> lambdaHelium;
    std::array<double, 31> lambdaNitrogen;
    std::array<double, 31> lambdaIron;

    // conex extensions
    int    resamplingMode;
    double modThreshold;
    double f19_cx;
    double f19_meson;
    double f19;

    bool m_has_extensions;

  public:

    // * standard fields fields

    int get_seed() const
    {return Seed1;}

    int get_particle() const
    {return Particle;}

    double get_alpha() const
    {return Alpha;}

    auto get_energy_min() const
    {return units::gigaelectron_volt_t<double>(util::math::pow(10, lgEmin - 9));}

    auto get_energy_max() const
    {return units::gigaelectron_volt_t<double>(util::math::pow(10, lgEmax - 9));}

    auto get_zenith_min() const
    {return units::degree_t<double>(zMin);}

    auto get_zenith_max() const
    {return units::degree_t<double>(zMax);}

    float get_version() const
    {return Version;}

    float get_output_version() const
    {return OutputVersion;}

    int get_he_model() const
    {return HEModel;}

    int get_le_model() const
    {return LEModel;}

    auto get_energy_threshold() const
    {return units::gigaelectron_volt_t<double>(HiLowEgy);}

    auto get_profile_cut_hadrons() const
    {return units::gigaelectron_volt_t<double>(hadCut);}

    auto get_profile_cut_em() const
    {return units::gigaelectron_volt_t<double>(emCut);}

    double get_threshold_hadrons() const
    {return hadThr;}

    double get_threshold_muons() const
    {return muThr;}

    double get_threshold_electrons() const
    {return emThr;}

    auto get_cut_hadrons() const
    {return units::gigaelectron_volt_t<double>(haCut);}

    auto get_cut_muons() const
    {return units::gigaelectron_volt_t<double>(muCut);}

    auto get_cut_electrons() const
    {return units::gigaelectron_volt_t<double>(elCut);}

    auto get_cut_photons() const
    {return units::gigaelectron_volt_t<double>(gaCut);}

    // mean free path arrays
    auto get_lambda_energy() const
    {
      using namespace units::literals;
      std::array<units::gigaelectron_volt_t<double>, 31> ret;
      for (int i = 0; i < 31; ++i) {
        ret[i] = std::pow(10, lambdaLgE[i] - 9) * 1_GeV;
      }
      return ret;
    }

    auto& get_lambda_proton() const
    {return reinterpret_cast<const std::array<units::grams_per_cubed_centimeter_t<double>, 31>&>(lambdaProton);}

    auto& get_lambda_pion() const
    {return reinterpret_cast<const std::array<units::grams_per_cubed_centimeter_t<double>, 31>&>(lambdaPion);}

    auto& get_lambda_helium() const
    {return reinterpret_cast<const std::array<units::grams_per_cubed_centimeter_t<double>, 31>&>(lambdaHelium);}

    auto& get_lambda_nitrogen() const
    {return reinterpret_cast<const std::array<units::grams_per_cubed_centimeter_t<double>, 31>&>(lambdaNitrogen);}

    auto& get_lambda_iron() const
    {return reinterpret_cast<const std::array<units::grams_per_cubed_centimeter_t<double>, 31>&>(lambdaIron);}

    // extension fields
    bool has_extensions() const
    {return m_has_extensions;}

    int get_resampling_mode() const
    {return resamplingMode;}

    double get_mod_threshold() const
    {
      // TODO
      std::cerr << "WARNING: modification threshold has unknown units and is computed as double" << std::endl;
      return modThreshold;
    }

    double get_f19_cx() const
    {return f19_cx;}

    double get_f19_meson() const
    {return f19_meson;}

    double get_f19() const
    {return f19;}
  };

} // namespace conex

#endif // _conex_header_h