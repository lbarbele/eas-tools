#ifndef _models_atmosphere_h
#define _models_atmosphere_h

#include <cmath>

namespace models::atmosphere {

  class us_standard {
  public:
  
    // * vertical mass overburden in g/cm^2 as a function of height in m
    double
    get_vertical_depth(
      const double height /* in m */
    ) const
    {
      if (height < 0) {
        return get_vertical_depth(0);
      } else if (height < 4e3) {
        return -186.5562 + 1222.6562 * std::exp(-height / 9941.8638);
      } else if (height < 1e4) {
        return -94.919 + 1144.9069 * std::exp(-height / 8781.5355);
      } else if (height < 4e4) {
        return 0.61289 + 1305.5948 * std::exp(-height / 6361.4304);
      } else if (height < 1e5) {
        return 0 + 540.1778 * std::exp(-height / 7721.7016);
      } else if (height < 1.128292e5) {
        return 0.01128292 - height / 1e7;
      } else {
        return 0;
      }
    }

    // * compute height im m given vertical mass overburden in g/cm^2 
    double
    get_height_from_depth(
      const double depth /* in g/cm^2 */
    ) const
    {
      if (depth <= 0) {
        return 1.128292e5;
      } else if (depth < 0.00128292) {
        return - 1e7 * (depth - 0.01128292);
      } else if (depth < 3.0395) {
        return - 7721.7016 * std::log(depth / 540.1778);
      } else if (depth < 271.7) {
        return - 6361.4304 * std::log((depth - 0.61289) / 1305.5948);
      } else if (depth < 631.101) {
        return - 8781.5355 * std::log((depth + 94.919) / 1144.9069);
      } else if (depth < 1036.1) {
        return - 9941.8638 * std::log((depth + 186.5562) / 1222.6562);
      } else {
        return get_height_from_depth(1036.1-0.0000001);
      }
      return 0;
    }

    // * air density in g/cm^2 as a function of height in m
    double
    get_density_from_height(
      const double height /* in m */
    ) const
    {
      if (height < 0) {
        return get_density_from_height(0);
      } else if (height < 4e3) {
        return (1222.6562 / 9941.8638) * std::exp(-height / 9941.8638);
      } else if (height < 1e4) {
        return (1144.9069 / 8781.5355) * std::exp(-height / 8781.5355);
      } else if (height < 4e4) {
        return (1305.5948 / 6361.4304) * std::exp(-height / 6361.4304);
      } else if (height < 1e5) {
        return ( 540.1778 / 7721.7016) * std::exp(-height / 7721.7016);
      } else if (height < 1.128292e5) {
        return 1.0 / 1e9;
      } else {
        return 0;
      }
    }

  };

}

#endif