#ifndef _util_units_h
#define _util_units_h

#include <units/units.h>
#include <units/angular.h>

namespace units {
  // * new units 

  // define grams per cm^2
  using grams_per_squared_centimeter = make_unit<gram, power<centimeter, -2>>;
  units_set_literal(grams_per_squared_centimeter, gcm2);
  units_set_quantity_alias(grams_per_squared_centimeter);

  // define grams per cm^3
  using grams_per_cubed_centimeter = make_unit<gram, power<centimeter, -3>>;
  units_set_literal(grams_per_cubed_centimeter, gcm3);
  units_set_quantity_alias(grams_per_cubed_centimeter);

  // * aliases for scalar quantities

  using depth_t = units::grams_per_squared_centimeter_t<double>;
  using density_t = units::grams_per_cubed_centimeter_t<double>;
  using length_t = units::meter_t<double>;
}

#endif
