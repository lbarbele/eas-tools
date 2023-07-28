#ifndef _models_atmosphere_h
#define _models_atmosphere_h

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string_view>
#include <ostream>

#include <util/math.h>
#include <util/point.h>
#include <util/vector.h>
#include <util/constants.h>

#include <util/units.h>

namespace models::atmosphere {
  
  using namespace units::literals;

  class us_standard {
  public:
    int m_nlayers;
    std::vector<units::depth_t> m_a;
    std::vector<units::depth_t> m_b;
    std::vector<units::length_t> m_c;
    std::vector<units::length_t> m_height_boundaries;
    std::vector<units::depth_t> m_depth_boundaries;

  public:

    us_standard()
    {
      m_nlayers = 5;

      m_a = {-186.555305_gcm2, -94.919_gcm2, 0.61289_gcm2, 0.0_gcm2, 0.01128292_gcm2};
      m_b = {1222.6562_gcm2, 1144.9069_gcm2, 1305.5948_gcm2, 540.1778_gcm2, 1_gcm2};
      m_c = {994186.38_cm, 878153.55_cm, 636143.04_cm, 772170.16_cm, 1e9_cm};

      m_height_boundaries = {0_m, 4e3_m, 1e4_m, 4e4_m, 1e5_m, 1.128292e5_m};
      m_depth_boundaries = {1.036101e3_gcm2, 6.311009e2_gcm2, 2.717009e2_gcm2, 3.039560_gcm2, 1.282922e-3_gcm2, 0_gcm2};
    }

    units::length_t max_height() const
    {return m_height_boundaries.back();}

    units::depth_t max_depth() const
    {return m_depth_boundaries.front();}

    // * get index of atmopsheric layer correspoding to given height
    int
    get_layer_index(
      units::length_t height
    ) const
    {
      if (height < m_height_boundaries.front()) {
        throw std::logic_error("height below boundary");
      }

      if (height > m_height_boundaries.back()) {
        throw std::logic_error("height above boundary");
      }

      for (int ilayer = 0; ilayer < m_nlayers; ++ilayer) {
        if (height <= m_height_boundaries[ilayer+1]) {
          return ilayer;
        }
      }

      // if here, there is some problem in m_height_boundaries
      throw std::logic_error("unable to find atmosphere layer from height");
    }

    // * same as above, but for given height
    int
    get_layer_index(
      const units::depth_t depth
    ) const
    {
      if (depth < m_depth_boundaries.back()) {
        throw std::logic_error("depth below boundary");
      }

      if (depth > m_depth_boundaries.front()) {
        throw std::logic_error("depth above boundary");
      }

      for (int ilayer = 0; ilayer < m_nlayers; ++ilayer) {
        if (depth >= m_depth_boundaries[ilayer+1]) {
          return ilayer;
        }
      }

      // if here, there is some problem in m_height_boundaries
      throw std::logic_error("unable to find atmosphere layer from depth");
    }
    
    // * vertical mass overburden as a function of height
    units::depth_t
    get_depth(
      const units::length_t h
    ) const
    {
      if (h > max_height()) {
        return units::depth_t(0);
      }

      const int ilayer = get_layer_index(h);
      return (ilayer < m_nlayers-1) ?
        m_a[ilayer] + m_b[ilayer] * std::exp(-h/m_c[ilayer]) :
        m_a[ilayer] - m_b[ilayer]*(h/m_c[ilayer]);
    }

    // * compute height given vertical mass overburden
    units::length_t
    get_height(
      const units::depth_t depth
    ) const
    {
      const int ilayer = get_layer_index(depth);
      return (ilayer < m_nlayers-1) ?
        m_c[ilayer] * std::log(m_b[ilayer]/(depth-m_a[ilayer])) :
        m_c[ilayer] * ((m_a[ilayer] - depth) / m_b[ilayer]);
    }

    // * air density as a function of height
    units::density_t
    get_density(
      const units::length_t h
    ) const
    {
      if (h > max_height()) {
        return units::density_t(0);
      }

      const int ilayer = get_layer_index(h);
      return (ilayer < m_nlayers-1)?
        (m_b[ilayer] / m_c[ilayer]) * std::exp(-h/m_c[ilayer]) : 
        (m_b[ilayer] / m_c[ilayer]);
    }

    // * air density as a function of depth
    units::density_t
    get_density(
      const units::depth_t depth
    ) const
    {
      const int ilayer = get_layer_index(depth);
      return (ilayer < m_nlayers-1)?
        (depth - m_a[ilayer])/m_c[ilayer] :
        m_b[ilayer]/m_c[ilayer];
    }

    // * air density at given point
    units::density_t
    get_density(
      const util::point_t<units::length_t>& p
    ) const
    {
      static const util::point_t<units::length_t> c(0_m, 0_m, -util::constants::earth_radius, util::frame::standard);
      const units::length_t height((p-c).norm() - util::constants::earth_radius);
      return get_density(height);
    }

    // * traversed mass between two points
    units::depth_t
    get_traversed_mass(
      const util::point_t<units::length_t>& a,
      const util::point_t<units::length_t>& b
    ) const
    {
      // point at earth's center
      static const util::point_t<units::length_t> earth_center(0_m, 0_m, -util::constants::earth_radius, util::frame::standard);

      // separation vector
      const util::vector_t<units::length_t> separation = b - a;

      // position vectors with origin at the earth's center
      const auto ra = a - earth_center;
      const auto rb = b - earth_center;

      // ensure integration path goes upwards
      const double cos_alpha = util::cos_angle(ra, separation);
      const double sin_alpha = std::sqrt((1+cos_alpha)*(1-cos_alpha));

      if (cos_alpha < 0) {
        // integration goes downwards, so check if there is a change of sign in cos_theta along the trajectory
        const units::length_t s_cross = - ra.norm() * cos_alpha;

        if (separation.norm() > s_cross) {
          // signal of cos_theta changes, so we split the integral in two parts
          const double f = s_cross / separation.norm();
          const util::point_t<units::length_t> cross_point_a = a + (1 - 1e-10) * f * separation;
          const util::point_t<units::length_t> cross_point_b = a + (1 + 1e-10) * f * separation;
          return get_traversed_mass(cross_point_a, a) + get_traversed_mass(cross_point_b, b);
        } else {
          // no change of signal in cos_theta, simply invert points
          return get_traversed_mass(b, a);
        }
      }

      // compute initial/final heights (in meters !)
      const units::length_t ha = ra.norm() - util::constants::earth_radius;
      const units::length_t hb = rb.norm() - util::constants::earth_radius;

      // get layer indices
      const auto ia = get_layer_index(ha);
      const auto ib = get_layer_index(hb);

      if (ia != ib) {
        // if integration path crosses layer boundaries, split the integral
        const units::length_t h_cross = m_height_boundaries[ia+1];
        const units::length_t r_cross = h_cross + util::constants::earth_radius;
        const units::length_t s_cross = util::math::sqrt((r_cross + ra.norm()*sin_alpha)*(r_cross - ra.norm()*sin_alpha)) - ra.norm()*cos_alpha;
        const double f = s_cross / separation.norm();
        const auto cross_point_a = a + (1 - 1e-10) * f * separation;
        const auto cross_point_b = a + (1 + 1e-10) * f * separation;
        return get_traversed_mass(a, cross_point_a) + get_traversed_mass(cross_point_b, b);
      }

      // here, we have a trajectory that goes upwards and does not cross any layer boundary
      const auto ilayer = ia;

      if (ilayer < m_nlayers-1) {
        const units::length_t catm = m_c[ilayer];

        const double lower = std::exp(-hb/catm);
        const double upper = std::exp(-ha/catm);
        const double tolerance = 1e-5;

        const units::length_t r0 = ra.cross_product(rb).norm() / separation.norm();

        const auto integrand = [=](const double u){
          const units::length_t r = util::constants::earth_radius - std::log(u) * catm;
          const double x = r0/r;
          return 1.0/std::sqrt((1+x)*(1-x));
        };

        return m_b[ilayer] * util::math::romberg_integral(lower, upper, tolerance, integrand);
      } else {
        return separation.norm() * m_b[ilayer] / m_c[ilayer];
      }
    }

    // * traversed length
    template <util::concepts::scalar U>
    util::point_t<units::length_t>
    propagate(
      const util::point_t<units::length_t>& initial_point,
      const util::vector_t<U>& direction_input,
      const units::depth_t traversed_mass
    ) const
    {
      // * constants, input conversions, and helpers
      
      // point at earth's center
      static const util::point_t<units::length_t> earth_center(0_m, 0_m, -util::constants::earth_radius, util::frame::standard);

      // parameters for the newton method below
      const units::length_t tolerance = 1_um; // 1 micrometer
      const unsigned max_iterations = 100;

      // normalize the direction and make it dimensionless
      const util::vector_d direction = direction_input.get_normalized(1);

      // * first approximation: assume X_slant = X_vertical / cos_theta

      // position vector (with origin at Earth's center)
      const util::vector_t<units::length_t> position_vector = initial_point - earth_center;

      // inclination of the trajectory (measured downwards!)
      const double cos_theta = -util::cos_angle(position_vector, direction);

      // initial altitude and depth
      const units::length_t starting_altitude = position_vector.norm() - util::constants::earth_radius;
      const units::depth_t starting_vertical_depth = get_depth(starting_altitude);

      // approximation to ending depth/altitude
      const units::depth_t ending_vertical_depth = starting_vertical_depth + traversed_mass * cos_theta;
      const units::length_t ending_altitude = get_height(ending_vertical_depth);

      // approximated displacement
      units::length_t displacement = (starting_altitude - ending_altitude) / cos_theta;

      // ending point
      util::point_t<units::length_t> ending_point = initial_point + displacement * direction;

      // * improve the solution using newton's method

      // loop up to max_iterations to find a valid solution
      for (unsigned i = 0; i < max_iterations; ++i) {
        // compute current step in the displacement (s), with change = (X(s) - X0) / X'(s), and X'(s) = rho(s)
        const units::length_t change = (get_traversed_mass(initial_point, ending_point) - traversed_mass) / get_density(ending_point);
        // compute new displacement
        displacement -= change;
        // compute new estimated ending point
        ending_point = initial_point + displacement * direction;
        // if change is below tolerance, we are done
        if (util::math::abs(change) < tolerance) {
          return ending_point;
        }
      }

      std::cerr << "unable to reach desired precision" << std::endl;
      throw;
    }

  };

}

#endif