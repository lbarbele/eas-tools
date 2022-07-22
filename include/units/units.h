#ifndef _units_units_h
#define _units_units_h

#include <units/units_impl.h>
#include <util/vector.h>

namespace units {
  // * define surface density as grams per square centimiter (same unit as atmospheric depth)
  UNIT_ADD(
    surface_density,
    gram_per_square_cm,
    grams_per_square_cm,
    gcm,
    compound_unit<mass::grams, inverse<squared<length::cm>>>
  )

  // * scalar quantities
  using density_t = density::kilograms_per_liter_t;
  using depth_t = surface_density::gram_per_square_cm_t;
  using height_t = length::meter_t;
  using time_t = time::second_t;

  // * vector quantities
  using position_t = util::vector<length::meter_t>;
}

using namespace units::literals;

#endif