#ifndef _util_constants_h
#define _util_constants_h

#include <util/units.h>

namespace util::constants {
  constexpr auto earth_radius = units::meter_t<double>(6371315);
  constexpr auto speed_of_light = units::speed_of_light_t<double>(1);
  constexpr auto pi = units::pi_t<double>(1);
  constexpr auto proton_mass = units::gev_per_c_squared_t<double>(0.93827208816);
} // namespace util::constants

#endif