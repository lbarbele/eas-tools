#ifndef _units_units_h
#define _units_units_h

#include <units/units_impl.h>

namespace units {
  // * define surface density as grams per square centimiter (same unit as atmospheric depth)
  UNIT_ADD(
    surface_density,
    gram_per_square_cm,
    grams_per_square_cm,
    gcm,
    compound_unit<mass::grams, inverse<units::squared<length::cm>>>
  )

  // * physical quantities using units
  using height_t = units::length::meter_t;
  using depth_t = units::surface_density::gram_per_square_cm_t;
  using density_t = units::density::kilograms_per_liter_t;
}

using namespace units::literals;

#endif