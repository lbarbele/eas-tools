#ifndef _models_atmosphere_h
#define _models_atmosphere_h

#include <cmath>
#include <vector>
#include <stdexcept>
#include <string_view>
#include <ostream>

#include <util/point.h>
#include <util/vector.h>
#include <util/constants.h>

namespace util::math {

	template<class Function>
	double
	romberg_integral (
    const double a,
    const double b,
    const double eps, 
    Function function
  ) {
		constexpr int max_it = 30;
		
		std::array<double, max_it> v = {0.0};
		
		double step = b-a;
		v[0] = 0.5*step*(function(a) + function(b));
		
		double value_before = v[0];
		int n_steps = 1;
		
		for (int k = 1; k < max_it; k++) {
			n_steps *= 2;
			step /= 2.0;
			
			for (int i = 1; i < n_steps; i+=2)
				v[k] += function(a + i*step);
			v[k] = v[k]*step + 0.5*v[k-1];
			
			for (int j = k-1; j >= 0; j--)
				v[j] = v[j+1] + (v[j+1] - v[j]) / (std::pow(4,k-j) - 1.0);
			
			if (std::fabs(value_before - v[0]) < eps)
				return v[0];
			else
				value_before = v[0];
		}
		
		return -1;
	}

}

namespace models::atmosphere {

  inline namespace units {

    template <class UnitT>
    class quantity_t {
    private:
      long double value;

    public:
      using type = quantity_t;
      using unit = UnitT;

      explicit constexpr quantity_t<UnitT>(const long double v = 0) : value(v) {}

      double get_value() const
      {return value;}


      constexpr type& operator *= (const long double x) {value *= x; return *this;}
      constexpr type& operator /= (const long double x) {value /= x; return *this;}
      constexpr type& operator += (const type other) {value += other.value; return *this;}
      constexpr type& operator -= (const type other) {value -= other.value; return *this;}

      constexpr type operator+(const type other) const {return type{value + other.value};}
      constexpr type operator-(const type other) const {return type{value - other.value};}

      constexpr type operator*(const long double x) const {return type{value*x};}
      constexpr type operator/(const long double x) const {return type{value/x};}

      constexpr long double operator/(const type other) const {return value/other.value;}

      constexpr type operator-() const {return type{-value};}
      constexpr type operator+() const {return type{+value};}

      constexpr bool operator<(const type other) const {return value < other.value;}
      constexpr bool operator>(const type other) const {return value > other.value;}

      constexpr bool operator<=(const type other) const {return value <= other.value;}
      constexpr bool operator>=(const type other) const {return value >= other.value;}

      constexpr bool operator==(const type other) const {return value == other.value;}
      constexpr bool operator!=(const type other) const {return value != other.value;}
    };

    // * multiply quantity from lhs
    template <class Un>
    constexpr auto operator*(const long double x, const quantity_t<Un>& q)
    {return q*x;}

    // * print quantity
    template <class CharT, class Traits, class UnitT>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os, const quantity_t<UnitT>& q)
    {return os << q.get_value() << " " << UnitT::symbol;}
    

    struct unit_t {
      static constexpr inline std::string_view symbol = "?";
    };

    struct gcm_t
    : public unit_t {
      static constexpr inline std::string_view symbol = "g cm^-2";
    };
    using depth_t = quantity_t<gcm_t>;

    struct meters_t
    : public unit_t {
      static constexpr inline std::string_view symbol = "m";
    };
    using height_t = quantity_t<meters_t>;

    struct gram_per_cubic_cm_t
    : public unit_t {
      static constexpr inline std::string_view symbol = "g cm^-3";
    };
    using density_t = quantity_t<gram_per_cubic_cm_t>;

  }

  inline namespace literals {
    constexpr auto operator"" _m(const long double value){return height_t{value};}
    constexpr auto operator"" _m(const unsigned long long value){return height_t{static_cast<long double>(value)};}

    constexpr auto operator"" _cm(const long double value){return height_t{0.01*value};}
    constexpr auto operator"" _cm(const unsigned long long value){return height_t{0.01*static_cast<long double>(value)};}

    constexpr auto operator"" _gcm2(const long double value){return depth_t{value};}
    constexpr auto operator"" _gcm2(const unsigned long long value){return depth_t{static_cast<long double>(value)};}

    constexpr auto operator"" _gcm3(const long double value){return density_t{value};}
    constexpr auto operator"" _gcm3(const unsigned long long value){return density_t{static_cast<long double>(value)};}
  }

  class us_standard {
  public:
    int m_nlayers;
    std::vector<depth_t> m_a;
    std::vector<depth_t> m_b;
    std::vector<height_t> m_c;
    std::vector<height_t> m_height_boundaries;
    std::vector<depth_t> m_depth_boundaries;

  public:

    us_standard()
    {
      m_nlayers = 5;

      m_a = {-186.555305_gcm2, -94.919_gcm2, 0.61289_gcm2, 0.0_gcm2, 0.01128292_gcm2};
      m_b = {1222.6562_gcm2, 1144.9069_gcm2, 1305.5948_gcm2, 540.1778_gcm2, 1_gcm2};
      m_c = {994186.38_cm, 878153.55_cm, 636143.04_cm, 772170.16_cm, 1e9_cm};

      m_height_boundaries = {0._m, 4e3_m, 1e4_m, 4e4_m, 1e5_m, 1.128292e5_m};

      m_depth_boundaries.resize(m_height_boundaries.size());
      for (int i = 0; i < m_nlayers; ++i) {
        m_depth_boundaries[i] = get_depth(m_height_boundaries[i]);
      }
    }

    // * get index of atmopsheric layer correspoding to given height
    int
    get_layer_index(
      height_t height
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
      const depth_t depth
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
    depth_t
    get_depth(
      const height_t h
    ) const
    {
      const int ilayer = get_layer_index(h);
      return (ilayer < m_nlayers-1) ?
        m_a[ilayer] + m_b[ilayer] * std::exp(-h/m_c[ilayer]) :
        m_a[ilayer] - m_b[ilayer]*(h/m_c[ilayer]);
    }

    // * compute height given vertical mass overburden
    height_t
    get_height(
      const depth_t depth
    ) const
    {
      const int ilayer = get_layer_index(depth);
      return (ilayer < m_nlayers-1) ?
        m_c[ilayer] * std::log(m_b[ilayer]/(depth-m_a[ilayer])) :
        m_c[ilayer] * ((m_a[ilayer] - depth) / m_b[ilayer]);
    }

    // * air density as a function of height
    density_t
    get_density(
      const height_t h
    ) const
    {
      const int ilayer = get_layer_index(h);
      return (ilayer < m_nlayers-1)?
        density_t((m_b[ilayer]/1_gcm2) / (m_c[ilayer]/1_cm)) * std::exp(-h/m_c[ilayer]) : 
        density_t((m_b[ilayer]/1_gcm2) / (m_c[ilayer]/1_cm));
    }

    // * air density as a function of depth
    density_t
    get_density(
      const depth_t depth
    ) const
    {
      const int ilayer = get_layer_index(depth);
      return (ilayer < m_nlayers-1)?
        density_t((depth/1_gcm2 - m_a[ilayer]/1_gcm2) / (m_c[ilayer]/1_cm)) : 
        density_t((m_b[ilayer]/1_gcm2) / (m_c[ilayer]/1_cm));
    }

    // * air density at given point
    density_t
    get_density(
      const util::point_t<height_t>& pi
    ) const
    {
      static const util::point_d c(0, 0, -util::constants::earth_radius, util::frame::standard);
      const util::point_d a(pi.x().get_value(), pi.y().get_value(), pi.z().get_value(), pi.get_frame());
      const height_t height((a-c).norm() - util::constants::earth_radius);
      return get_density(height);
    }

    // * traversed mass between two points
    depth_t
    get_traversed_mass(
      const util::point_t<height_t>& ai,
      const util::point_t<height_t>& bi
    ) const
    {
      // alias to earth radius
      static const double rea = util::constants::earth_radius;

      // point at earth's center
      static const util::point_d c(0, 0, -rea, util::frame::standard);

      // convert input to double (in meters!!!)
      const util::point_d a(ai.x().get_value(), ai.y().get_value(), ai.z().get_value(), ai.get_frame());
      const util::point_d b(bi.x().get_value(), bi.y().get_value(), bi.z().get_value(), bi.get_frame());

      // separation vector
      const util::vector_d separation = b - a;

      // position vectors with origin at the earth's center
      const auto ra = a - c;
      const auto rb = b - c;

      // check if integration path goes upwards
      double cos_alpha = ra.get_normalized(1) * separation.get_normalized(1);
      double sin_alpha = std::sqrt((1+cos_alpha)*(1-cos_alpha));

      if (cos_alpha < 0) {
        // integration goes downwards, so check if there is a change of sign in cos_theta along the trajectory
        const double s_cross = - ra.norm() * cos_alpha;

        if (separation.norm() > s_cross) {
          // signal of cos_theta changes, so we split the integral in two parts
          const double f = s_cross / separation.norm();
          const auto cpt_a = a + (1 - 1e-10) * f * separation;
          const auto cpt_b = a + (1 + 1e-10) * f * separation;
          const util::point_t<height_t> cross_point_a(cpt_a.x() * 1_m, cpt_a.y() * 1_m, cpt_a.z() * 1_m, cpt_a.get_frame());
          const util::point_t<height_t> cross_point_b(cpt_b.x() * 1_m, cpt_b.y() * 1_m, cpt_b.z() * 1_m, cpt_b.get_frame());
          return get_traversed_mass(cross_point_a, ai) + get_traversed_mass(cross_point_b, bi);
        } else {
          // no change of signal in cos_theta, simply invert points
          return get_traversed_mass(bi, ai);
        }
      }

      // compute initial/final heights (in meters !)
      const height_t ha(ra.norm() - rea);
      const height_t hb(rb.norm() - rea);

      // get layer indices
      const auto ia = get_layer_index(ha);
      const auto ib = get_layer_index(hb);

      if (ia != ib) {
        // if integration path crosses layer boundaries, split the integral
        const double h_cross = m_height_boundaries[ia+1] / 1_m;
        const double r_cross = h_cross + rea;
        const double s_cross = std::sqrt((r_cross + ra.norm()*sin_alpha)*(r_cross - ra.norm()*sin_alpha)) - ra.norm()*cos_alpha;
        const double f = s_cross / separation.norm();
        const auto cpt_a = a + (1 - 1e-10) * f * separation;
        const auto cpt_b = a + (1 + 1e-10) * f * separation;
        const util::point_t<height_t> cross_point_a(cpt_a.x() * 1_m, cpt_a.y() * 1_m, cpt_a.z() * 1_m, cpt_a.get_frame());
        const util::point_t<height_t> cross_point_b(cpt_b.x() * 1_m, cpt_b.y() * 1_m, cpt_b.z() * 1_m, cpt_b.get_frame());
        return get_traversed_mass(ai, cross_point_a) + get_traversed_mass(cross_point_b, bi);
      }

      // here, we have a trajectory that goes upwards and does not cross any layer boundary
      const auto ilayer = ia;

      if (ilayer < m_nlayers-1) {
        const double lower = std::exp(-hb/m_c[ilayer]);
        const double upper = std::exp(-ha/m_c[ilayer]);
        const double tolerance = 1e-8;

        const double catm = (m_c[ilayer] / 1_m);
        const double r0 = ra.cross_product(rb).norm() / separation.norm();

        const auto integrand = [=](const double u){
          const double r = rea - std::log(u) * catm;
          const double x = r0/r;
          return 1.0/std::sqrt((1+x)*(1-x));
        };

        double integral = util::math::romberg_integral(lower, upper, tolerance, integrand);
        integral *= m_b[ilayer] / 1_gcm2;
      
        return depth_t(integral);
      } else {
        return depth_t(separation.norm() * 100 * (get_density(ha) / 1_gcm3));
      }
    }

    // * traversed length
    util::point_t<height_t>
    propagate(
      const util::point_t<height_t>& ai,
      const util::vector_t<height_t>& d,
      const depth_t traversed_mass
    ) const
    {
      // - constants and helpers
      
      static const double rea = util::constants::earth_radius;
      static const util::point_d earth_center(0, 0, -rea, util::frame::standard);

      const double tolerance = 1e-6;
      const unsigned max_iterations = 100;

      const auto to_meter = [](const util::point_d p) {
        return util::point_t<height_t>(p.x()*1_m, p.y()*1_m, p.z()*1_m, p.get_frame());
      };

      // - unit conversion

      const auto initial_point = util::point_d(ai.x()/1_m, ai.y()/1_m, ai.z()/1_m, ai.get_frame());
      const auto direction = util::vector_d(d.x()/1_m, d.y()/1_m, d.z()/1_m, d.get_frame()).get_normalized(1);

      // - first approximation: assume X_slant = X_vertical / cos_theta

      // * position vector
      const util::vector_d position_vector = initial_point - earth_center;

      // * inclination of the trajectory (measured downwards!)
      const double cos_theta = -util::cos_angle(position_vector, direction);

      if (std::fabs(cos_theta) < 1e-5) {
        std::cerr << "bad cos_theta" << std::endl;
        throw;
      }

      // * initial altitude and depth
      const height_t starting_altitude = 1_m * (position_vector.norm() - rea);
      const depth_t starting_vertical_depth = get_depth(starting_altitude);

      // * approximation to ending depth/altitude
      const depth_t ending_vertical_depth = starting_vertical_depth + traversed_mass * cos_theta;
      const height_t ending_altitude = get_height(ending_vertical_depth);

      // * approximated displacement
      double displacement = ((starting_altitude - ending_altitude) / 1_m ) / cos_theta;

      // * ending point
      auto estimated_ending_point = initial_point + displacement * direction;

      // - improve the solution using newton's method

      // * loop up to max_iterations to find a valid solution
      for (unsigned i = 0; i < max_iterations; ++i) {
        // convert current estimate to a point carrying units
        const auto bi = to_meter(estimated_ending_point);
        // compute current step in the displacement (s), with change = (X(s) - X0) / X'(s), and X'(s) = rho(s)
        const double change = 0.01 * ((get_traversed_mass(ai, bi)-traversed_mass)/1_gcm2) / (get_density(bi)/1_gcm3);
        // compute new displacement
        displacement -= change;
        // compute new estimated ending point
        estimated_ending_point = initial_point + displacement * direction;
        // if change is below tolerance, we are done
        if (std::fabs(change) < tolerance) {
          return to_meter(estimated_ending_point);
        }
      }

      std::cerr << "unable to reach desired precision" << std::endl;
      throw;
    }

  };

}

#endif