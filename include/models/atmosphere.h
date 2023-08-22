#ifndef _models_atmosphere_h
#define _models_atmosphere_h

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <util/math.h>
#include <util/point.h>
#include <util/vector.h>
#include <util/constants.h>

#include <util/units.h>

namespace models::atmosphere {
  
  using namespace units::literals;

  class us_standard {
  private:

    struct layer {
      enum type {
        linear,
        exponential
      };

      const uint id;
      const units::depth_t a;
      const units::depth_t b;
      const units::length_t c;
      const units::length_t min_height;
      const units::length_t max_height;
      const units::depth_t min_depth;
      const units::depth_t max_depth;
      const type shape;
    };

    const std::array<layer, 5> m_layers = {
      //    id,               a,              b,           c,  hlow,          hup,            xlow,             xup, type
      layer{1u,  -186.5562_gcm2, 1222.6562_gcm2, 9941.8638_m,   0_m,        4e3_m,      631.1_gcm2,   1036.101_gcm2, layer::exponential},
      layer{2u,    -94.919_gcm2, 1144.9069_gcm2, 8781.5355_m, 4e3_m,        1e4_m,      271.7_gcm2,      631.1_gcm2, layer::exponential},
      layer{3u,    0.61289_gcm2, 1305.5948_gcm2, 6361.4304_m, 1e4_m,        4e4_m,     3.0396_gcm2,      271.7_gcm2, layer::exponential},
      layer{4u,        0.0_gcm2,  540.1778_gcm2, 7721.7016_m, 4e4_m,        1e5_m, 0.00128292_gcm2,     3.0396_gcm2, layer::exponential},
      layer{5u, 0.01128292_gcm2,       1.0_gcm2,       1e7_m, 1e5_m, 1.128292e5_m,          0_gcm2, 0.00128292_gcm2, layer::linear}
    };

    // * implementation of get_traversed mass

    units::depth_t
    do_get_traversed_mass(
      const units::length_t xa,
      const units::length_t xb,
      const units::length_t y
    ) const
    {
      const auto rea = util::constants::earth_radius;

      // we assume here the path mapped by xa and xb is never passing below ground!
      // this condition must (and, in fact, it is) ensured in get_traversed_mass()
      const auto height = [=](const auto x){
        return util::math::max(util::math::hypot(x, y) - rea, 0_m);
      };

      // the returning value
      units::depth_t traversed_mass(0.);

      // loop over layers
      const auto ia = get_layer(height(xa)).id;
      const auto ib = get_layer(height(xb)).id;

      const auto integrand = [=](const double u, const units::length_t c){
        const auto st = y/(rea - c*std::log(u));
        const auto ct = util::math::sqrt((1+st)*(1-st));
        return 1.0/ct;
      };

      for (uint idlay = ia; idlay <= ib; ++idlay) {
        const auto l = get_layer(idlay);

        const auto rl = rea + l.min_height;
        const auto ru = rea + l.max_height;

        const auto xl = idlay == ia ? xa : util::math::sqrt((rl+y)*(rl-y));
        const auto xu = idlay == ib ? xb : util::math::sqrt((ru+y)*(ru-y));

        if (l.shape == layer::linear) {
          // linear evolution
          traversed_mass += (xu - xl)*(l.b/l.c);
        } else if ((xu - xl)/l.c < 1e-7) {
          // exponential layer, but small displacement: linear evolution
          // the next term is proportional to ((xu-xl)/c)^2
          // so the error in this approximation is < ((xu-xl)/c)^2
          const auto xm = 0.5*(xl + xu); // middle point
          const auto hm = util::math::hypot(xm, y) - rea; // height at the middle point
          traversed_mass = get_density(hm) * (xu - xl); // linear approximation
        } else {
          const auto lower = util::math::exp(-height(xu)/l.c);
          const auto upper = util::math::exp(-height(xl)/l.c);

          try {
            traversed_mass += l.b * util::math::romberg_integral(lower, upper, 1e-8, integrand, l.c);
          } catch (...) {
            std::cerr << "failed to compute traversed mass" << std::endl;
            std::cerr << "layer id: " << idlay << std::endl;
            std::cerr << "xl:       " << xl << std::endl;
            std::cerr << "xu:       " << xu << std::endl;
            std::cerr << "y:        " << y << std::endl;
            throw std::runtime_error("this is a wrapped exception");
          }
        }
      }

      return traversed_mass;
    }

    // * determine the atmospheric layer corresponding to height, depth, or position

    const layer&
    get_layer(
      const unsigned int id
    ) const
    {
      for (const auto& layer : m_layers) {
        if (layer.id == id) {
          return layer;
        }
      }

      std::cerr
        << "unable to find atmosphere layer from id" << std::endl
        << "id was " << id << std::endl;

      throw std::runtime_error("atmosphere error");
    }

    const layer&
    get_layer(
      const units::length_t height
    ) const
    {
      if (height == max_height()) {
        return m_layers.back();
      }

      for (const auto& layer : m_layers) {
        if (layer.min_height <= height && height < layer.max_height) {
          return layer;
        }
      }

      std::cerr
        << "unable to find atmosphere layer from height" << std::endl
        << "height was " << height.convert<units::meter>() << std::endl;

      throw std::runtime_error("atmosphere error");
    }

    const layer&
    get_layer(
      const units::depth_t depth
    ) const
    {
      if (depth == min_depth()) {
        return m_layers.back();
      }

      for (const auto& layer : m_layers) {
        if (layer.min_depth < depth && depth <= layer.max_depth) {
          return layer;
        }
      }

      std::cerr
        << "unable to find atmosphere layer from depth" << std::endl
        << "depth was " << depth/1_gcm2 << " g/cm²" << std::endl;

      throw std::runtime_error("atmosphere error");
    }

    const layer& get_layer(const util::point_t<units::length_t>& position) const
    {return get_layer(get_height(position));}

  public:

    // * atmosphere boundaries

    units::length_t max_height() const
    {return m_layers.back().max_height;}

    units::length_t min_height() const
    {return m_layers.front().min_height;}

    units::depth_t min_depth() const
    {return m_layers.back().min_depth;}

    units::depth_t max_depth() const
    {return m_layers.front().max_depth;}
    
    // * compute vertical mass overburden at position

    units::depth_t
    get_depth(
      const units::length_t height
    ) const
    {
      if (height >= max_height()) {
        return units::depth_t(0);
      }

      const auto& l = get_layer(height);

      return l.shape == layer::exponential ? 
        l.a + l.b * util::math::exp(-height/l.c) :
        l.a - l.b * height / l.c;
    }

    units::depth_t get_depth(const util::point_t<units::length_t>& position) const
    {return get_depth(get_height(position));}

    // * compute height at given position

    units::length_t
    get_height(
      const units::depth_t depth
    ) const
    {
      const auto& l = get_layer(depth);

      const auto h = l.shape == layer::exponential ? 
        l.c * std::log(l.b / (depth - l.a)) :
        l.c * (l.a - depth) / l.b;

      return util::math::max(0_m, h);
    }

    units::length_t
    get_height(
      const util::point_t<units::length_t>& position 
    ) const
    {
      const auto earth_radius = util::constants::earth_radius;
      const auto earth_center = util::point_t<units::length_t>(0_m, 0_m, -earth_radius, util::frame::standard);
      return (position - earth_center).norm() - earth_radius;
    }

    // * compute air density at the given position

    units::density_t
    get_density(
      const units::length_t height
    ) const
    {
      if (height >= max_height()) {
        return units::density_t(0.);
      }

      const auto& l = get_layer(height);

      return l.shape == layer::exponential ?
        (l.b/l.c) * util::math::exp(-height/l.c) :
        (l.b/l.c);
    }

    units::density_t
    get_density(
      const units::depth_t depth
    ) const
    {
      const auto& l = get_layer(depth);

      return l.shape == layer::exponential ?
        (depth - l.a)/l.c :
        l.b/l.c;
    }

    units::density_t get_density(const util::point_t<units::length_t>& pos) const
    {return get_density(get_height(pos));}

    // * traversed mass between two points

    units::depth_t
    get_traversed_mass(
      const util::point_t<units::length_t>& a,
      const util::point_t<units::length_t>& b
    ) const
    {
      const auto rea = util::constants::earth_radius;
      const auto rmx = rea + max_height() - 1_mm;
      const auto cnt = util::point_t<units::length_t>(0_m, 0_m, -rea, util::frame::standard);

      // separation vector (norm and direction)
      const auto dir = (b - a).get_normalized(1.);

      // position vectors wrt to earth center
      const auto ra = a - cnt;
      const auto rb = b - cnt;

      if (ra.norm()-rea < -1_nm || rb.norm()-rea < -1_nm) {
        std::cerr
          << "error: tried to compute traversed mass for a point below ground level\n"
          << "ha: " << ra.norm() - rea << '\n'
          << "hb: " << rb.norm() - rea << '\n';
    
        throw std::runtime_error("atmosphere: get_traversed_mass error");
      }

      // impact radius
      const auto y = dir.cross_product(ra).norm();

      // positions along the integration direction
      const auto xm = util::math::sqrt((rmx+y)*(rmx-y));
      const auto xa = util::math::max(ra*dir, -xm);
      const auto xb = util::math::min(rb*dir, +xm);

      if (xa >= xm || xb <= -xm) {
        return units::depth_t(0.);
      } else if (xb < util::math::max(0_m, -xa)) {
        return do_get_traversed_mass(-xb, -xa, y);
      } else if (xa < 0_m) {
        return 2*do_get_traversed_mass(0_m, -xa, y) + do_get_traversed_mass(-xa, xb, y);
      } else {
        return do_get_traversed_mass(xa, xb, y);
      }
    }

    // * compute slant depth along z axis

    template <util::concepts::scalar U>
    units::depth_t
    get_slant_depth(
      const util::vector_t<U>& input_axis,
      const units::length_t distance,
      const units::length_t ground_level = 0_m
    ) const
    {
      // the integration axis and its sine/cosine directions
      const auto axis = input_axis.on_frame(util::frame::standard).get_normalized(1.);
      const auto cosa = axis.z();
      const auto sina = util::math::sqrt((1+cosa)*(1-cosa));

      // the distance from ground to atmosphere boundary along the given axis
      const auto rmax = util::constants::earth_radius + max_height();
      const auto rgnd = util::constants::earth_radius + ground_level;
      const auto dmax = util::math::sqrt((rmax+rgnd*sina)*(rmax-rgnd*sina)) - rgnd*cosa;

      // initial and final positions for integration
      const auto grnd = util::point_t<units::length_t>(0_m, 0_m, ground_level, util::frame::standard);
      const auto pini = grnd + axis*dmax;
      const auto pend = grnd + axis*distance;

      return get_traversed_mass(pini, pend);
    }

    // * compute traversed length in given direction

    template <util::concepts::scalar U>
    util::point_t<units::length_t>
    transport(
      const util::point_t<units::length_t>& initial_position,
      const util::vector_t<U>& input_direction,
      const units::depth_t& depth_interval,
      const units::length_t ground_level = 0_m,
      const units::length_t tolerance = 1_nm
    ) const
    {
      const auto rea = util::constants::earth_radius;
      const auto earth_center = util::point_t<units::length_t>(0_m, 0_m, -rea, util::frame::standard);

      // atmoisphere  limits considering the given ground level
      const auto hmin = ground_level;
      const auto hmax = max_height();

      const auto zmin = min_depth();
      const auto zmax = get_depth(hmin);

      auto direction = input_direction.get_normalized(1.);
      auto position = initial_position;
      auto dz = depth_interval;

      for (uint iter = 0; iter < 100; ++iter) {
        const auto cosine = util::cos_angle(position - earth_center, direction);

        // * compute approximation to displacement corresponding to dz

        units::length_t displacement = 0_m;

        if (util::math::abs(cosine) < 1e-2) {
          // case 1: (almost) horizontal displacement
          displacement = dz/get_density(position);
        } else {
          // case 2: non-horizontal displacement
          const auto final_depth = get_depth(position) - dz*cosine;

          const auto initial_height = get_height(position);

          const auto final_height =
            final_depth <= zmin ? hmax :
            final_depth >= zmax ? hmin :
            get_height(final_depth);

          displacement = (final_height - initial_height) / cosine;
        }

        // * update position and depth interval

        // update position
        const auto previous_position = position;
        position += displacement*direction;

        // don't allow heights below the ground level
        // (may happen due to floating point (im)precision)
        while (get_height(position) < ground_level) {
          position -= tolerance*direction;
        }

        // if displacement is smaller than tolerance, we are done
        if (util::math::abs(displacement) < tolerance) {
          return position;
        }

        // update depth
        dz -= get_traversed_mass(position, previous_position);

        // never use negative dz, change direction sign instead
        if (dz < 0_gcm2) {
          dz = -dz;
          direction = -direction;
        }
      }

      throw std::runtime_error("atm transport was unable to reach the desired precision");
    }

  };

}

#endif