#ifndef _conex_header_h
#define _conex_header_h

#include <array>

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
    // basic fields
    int get_seed() const
    {return Seed1;}

    int get_particle() const
    {return Particle;}

    double get_alpha() const
    {return Alpha;}

    double get_lge_min() const
    {return lgEmin;}

    double get_lge_max() const
    {return lgEmax;}

    double get_zenith_min() const
    {return zMin;}

    double get_zenith_max() const
    {return zMax;}

    float get_version() const
    {return Version;}

    float get_output_version() const
    {return OutputVersion;}

    int get_he_model() const
    {return HEModel;}

    int get_le_model() const
    {return LEModel;}

    float get_energy_threshold() const
    {return HiLowEgy;}

    float get_profile_cut_hadrons() const
    {return hadCut;}

    float get_profile_cut_em() const
    {return emCut;}

    float get_threshold_haddrons() const
    {return hadThr;}

    float get_threshold_muons() const
    {return muThr;}

    float get_threshold_electrons() const
    {return emThr;}

    float get_cut_hadrons() const
    {return haCut;}

    float get_cut_muons() const
    {return muCut;}

    float get_cut_electrons() const
    {return elCut;}

    float get_cut_photons() const
    {return gaCut;}

    // mean free path arrays
    const std::array<double, 31>& get_lambda_lge() const 
    {return lambdaLgE;}

    const std::array<double, 31>& get_lambda_proton() const
    {return lambdaProton;}

    const std::array<double, 31>& get_lambda_pion() const
    {return lambdaPion;}

    const std::array<double, 31>& get_lambda_helium() const
    {return lambdaHelium;}

    const std::array<double, 31>& get_lambda_nitrogen() const
    {return lambdaNitrogen;}

    const std::array<double, 31>& get_lambda_iron() const
    {return lambdaIron;}

    // extension fields
    bool has_extensions() const
    {return m_has_extensions;}

    int get_resampling_mode() const
    {return resamplingMode;}

    double get_mod_threshold() const
    {return modThreshold;}

    double get_f19_cx() const
    {return f19_cx;}

    double get_f19_meson() const
    {return f19_meson;}

    double get_f19() const
    {return f19;}
  };

} // namespace conex

#endif // _conex_header_h