#ifndef _util_constants_h
#define _util_constants_h

#include <util/units.h>

namespace util::constants {
  // constants in SI
  constexpr double c = 299792458.; // m/s
  constexpr double pi = 3.141592653589793238462643;
  constexpr auto earth_radius = units::meter_t<double>(6371315);

  // other units
  constexpr double proton_mass_gev = 0.93827208816;

} // namespace util::constants

#endif