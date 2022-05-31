#include <models/dedx_profile.h>

namespace models::dedx_profile {

  double
  get_usp_l(
    const double lgecal,
    const unsigned int id
  )
  {
    switch (id) {
      case  100: return 227.183 + lgecal*(7.16650 + 0.950551*lgecal);
      case  400: return 229.021 + lgecal*(5.90094 + 0.523577*lgecal);
      case 1200: return 228.763 + lgecal*(5.09137 + 0.366554*lgecal);
      case 2800: return 228.355 + lgecal*(4.70160 + 0.220192*lgecal);
      case 5600: return 227.620 + lgecal*(4.50438 + 0.175343*lgecal);
    }
    return 0;
  }

  double
  get_usp_r(
    const double lgecal,
    const unsigned int id
  )
  {
    switch (id) {
      case  100: return 0.256203 - lgecal*(0.0299802 - 0.00379108*lgecal);
      case  400: return 0.272744 - lgecal*(0.0329877 - 0.00371185*lgecal);
      case 1200: return 0.290858 - lgecal*(0.0362762 - 0.00336855*lgecal);
      case 2800: return 0.305702 - lgecal*(0.0390715 - 0.00327661*lgecal);
      case 5600: return 0.317921 - lgecal*(0.0411374 - 0.00328272*lgecal);
    }
    return 0;
  }

  double
  get_gh_x0(
    const double lgecal,
    const unsigned int id
  )
  {
    switch (id) {
      case  100: return -122.225 - lgecal*(69.6578 + 4.57977*lgecal);
      case  400: return -112.308 - lgecal*(60.9803 + 4.84916*lgecal);
      case 1200: return -89.7777 - lgecal*(53.2338 + 7.11282*lgecal);
      case 2800: return -72.9420 - lgecal*(49.0044 + 7.38145*lgecal);
      case 5600: return -60.2313 - lgecal*(44.8201 + 7.62627*lgecal);
    }
    return 0;
  }

  double
  get_gh_lambda(
    const double lgecal,
    const unsigned int id
  )
  {
    switch (id) {
      case  100: return 58.7702 - lgecal*(4.95382 - 0.880473*lgecal);
      case  400: return 62.7939 - lgecal*(5.91256 - 0.770987*lgecal);
      case 1200: return 66.6921 - lgecal*(6.70239 - 0.625783*lgecal);
      case 2800: return 69.8763 - lgecal*(7.38438 - 0.590561*lgecal);
      case 5600: return 72.4171 - lgecal*(7.85206 - 0.581305*lgecal);
    }
    return 0;
  }

} // namespace models::dedx_profile