#ifndef _util_units_h
#define _util_units_h

#include <units/angular.h>
#include <units/framework.h>
#include <units/units.h>

namespace units {
  // * new units 

  // define meter per second
  using meter_per_second = make_unit<meter, inverse<second>>;
  units_set_quantity_alias(meter_per_second);

  // define grams per cm^2
  using grams_per_squared_centimeter = make_unit<gram, power<centimeter, -2>>;
  units_set_literal(grams_per_squared_centimeter, gcm2);
  units_set_quantity_alias(grams_per_squared_centimeter);

  // define grams per cm^3
  using grams_per_cubed_centimeter = make_unit<gram, power<centimeter, -3>>;
  units_set_literal(grams_per_cubed_centimeter, gcm3);
  units_set_quantity_alias(grams_per_cubed_centimeter);

  // energy in electron volt
  units_add_derived_unit(electron_volt, eV, make_unit<joule, yocto, ratio<1'602'176'634, 10'000>>);
  units_set_prefixes(electron_volt, eV, kilo, mega, giga, tera, peta, exa, zetta);

  // speed of light
  units_add_derived_unit(speed_of_light, c, make_unit<ratio<299'792'458>, meter, inverse<second>>);

  // gev per c
  using gev_per_c = make_unit<gigaelectron_volt, inverse<speed_of_light>>;
  units_set_quantity_alias(gev_per_c);
  units_set_symbol(gev_per_c, GeV c^-1);

  // gev per c^2
  using gev_per_c_squared = make_unit<gev_per_c, inverse<speed_of_light>>;
  units_set_quantity_alias(gev_per_c_squared);
  units_set_symbol(gev_per_c_squared, GeV c^-2);

  // inverse grams per squared centimeter
  using squared_centimeters_per_gram = make_unit<inverse<grams_per_squared_centimeter>>;
  units_set_quantity_alias(squared_centimeters_per_gram);

  // giga electron volt per (g cm^-2)
  using gev_per_gcm = make_unit<gigaelectron_volt, inverse<grams_per_squared_centimeter>>;
  units_set_quantity_alias(gev_per_gcm);

  // * aliases for scalar quantities

  using depth_t = grams_per_squared_centimeter_t<double>;
  using density_t = grams_per_cubed_centimeter_t<double>;
  using length_t = meter_t<double>;
  using energy_t = gigaelectron_volt_t<double>;
  using momentum_t = quantity<make_unit<gigaelectron_volt, inverse<speed_of_light>>, double>;
  using time_t = second_t<double>;
  using speed_t = quantity<make_unit<meter, inverse<second>>, double>;
  using angle_t = radian_t<double>;
  using mass_t = quantity<make_unit<gigaelectron_volt, inverse_squared<speed_of_light>>, double>;
  using energy_deposit_t = gev_per_gcm_t<double>;
}

#endif
