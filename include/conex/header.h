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
    int                                Seed1;
    int                                Particle;
    double                             Alpha;
    double                             lgEmin;
    double                             lgEmax;
    units::degree_t<double>            zMin;
    units::degree_t<double>            zMax;
    float                              Version;
    float                              OutputVersion;
    int                                HEModel;
    int                                LEModel;
    units::gigaelectron_volt_t<float>  HiLowEgy;
    units::gigaelectron_volt_t<float>  hadCut;
    units::gigaelectron_volt_t<float>  emCut;
    float                              hadThr;
    float                              muThr;
    float                              emThr;
    units::gigaelectron_volt_t<float>  haCut;
    units::gigaelectron_volt_t<float>  muCut;
    units::gigaelectron_volt_t<float>  elCut;
    units::gigaelectron_volt_t<float>  gaCut;
    
    std::array<double, 31> lambdaLgE;
    std::array<units::grams_per_squared_centimeter_t<double>, 31> lambdaProton;
    std::array<units::grams_per_squared_centimeter_t<double>, 31> lambdaPion;
    std::array<units::grams_per_squared_centimeter_t<double>, 31> lambdaHelium;
    std::array<units::grams_per_squared_centimeter_t<double>, 31> lambdaNitrogen;
    std::array<units::grams_per_squared_centimeter_t<double>, 31> lambdaIron;

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

    units::energy_t get_energy_min() const
    {return units::gigaelectron_volt_t<double>(util::math::pow(10, lgEmin - 9));}

    units::energy_t get_energy_max() const
    {return units::gigaelectron_volt_t<double>(util::math::pow(10, lgEmax - 9));}

    units::angle_t get_zenith_min() const
    {return zMin;}

    units::angle_t get_zenith_max() const
    {return zMax;}

    float get_version() const
    {return Version;}

    float get_output_version() const
    {return OutputVersion;}

    int get_he_model() const
    {return HEModel;}

    int get_le_model() const
    {return LEModel;}

    units::energy_t get_energy_threshold() const
    {return HiLowEgy;}

    units::energy_t get_profile_cut_hadrons() const
    {return hadCut;}

    units::energy_t get_profile_cut_em() const
    {return emCut;}

    double get_threshold_hadrons() const
    {return hadThr;}

    double get_threshold_muons() const
    {return muThr;}

    double get_threshold_electrons() const
    {return emThr;}

    units::energy_t get_cut_hadrons() const
    {return haCut;}

    units::energy_t get_cut_muons() const
    {return muCut;}

    units::energy_t get_cut_electrons() const
    {return elCut;}

    units::energy_t get_cut_photons() const
    {return gaCut;}

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

    const auto& get_lambda_proton() const
    {return lambdaProton;}

    const auto& get_lambda_pion() const
    {return lambdaPion;}

    const auto& get_lambda_helium() const
    {return lambdaHelium;}

    const auto& get_lambda_nitrogen() const
    {return lambdaNitrogen;}

    const auto& get_lambda_iron() const
    {return lambdaIron;}

    // extension fields
    bool has_extensions() const
    {return m_has_extensions;}

    int get_resampling_mode() const
    {return resamplingMode;}

    units::energy_t get_mod_threshold() const
    {return units::gigaelectron_volt_t<double>(std::pow(10, modThreshold));}

    double get_f19_cx() const
    {return f19_cx;}

    double get_f19_meson() const
    {return f19_meson;}

    double get_f19() const
    {return f19;}
  };

} // namespace conex

#endif // _conex_header_h