#ifndef _models_dedx_profile_h
#define _models_dedx_profile_h

namespace models::dedx_profile {

  // parametrizations obtained using the tool fit_conex_profiles (see the 
  // cpp file for details) using 120000 showers for each primary simulated with 
  // sibyll 2.3d. The energy range of the simulations is 10^17 to 10^20 eV.
  double get_usp_l(const double lgecal, const unsigned int conex_id);
  double get_usp_r(const double lgecal, const unsigned int conex_id);
  double get_gh_x0(const double lgecal, const unsigned int conex_id);
  double get_gh_lambda(const double lgecal, const unsigned int conex_id);

} // namespace models::dedx_profile

#endif 