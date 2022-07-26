#ifndef _models_atmosphere_h
#define _models_atmosphere_h

#include <cmath>
#include <vector>
#include <stdexcept>

#include <util/vector.h>
#include <util/constants.h>

namespace models::atmosphere {

  template <class UnitT>
  class quantity_t {
  private:
    long double value;

  public:
    using type = quantity_t<UnitT>;

    constexpr quantity_t<UnitT>() {}
    constexpr quantity_t<UnitT>(const long double v) : value(v) {}


    inline constexpr type& operator *= (const long double x) {value *= x; return this;}
    inline constexpr type& operator /= (const long double x) {value *= x; return this;}


    inline constexpr type operator+(const type other) const {return type{value + other.value};}
    inline constexpr type operator-(const type other) const {return type{value - other.value};}

    inline constexpr type operator*(const long double x) const {return type{x*value};}
    inline constexpr type operator/(const long double x) const {return type{x/value};}

    inline constexpr long double operator/(const type other) const {return value/other.value;}

    inline constexpr type operator-() const {return type{-value};}
    inline constexpr type operator+() const {return type{+value};}

    inline constexpr bool operator<(const type other) const {return value < other.value;}
    inline constexpr bool operator>(const type other) const {return value > other.value;}

    inline constexpr bool operator<=(const type other) const {return value <= other.value;}
    inline constexpr bool operator>=(const type other) const {return value >= other.value;}

    inline constexpr bool operator==(const type other) const {return value == other.value;}
    inline constexpr bool operator!=(const type other) const {return value != other.value;}
  };

  struct unit_t {};

  struct gcm_t : public unit_t {};
  using depth_t = quantity_t<gcm_t>;

  struct meters_t : public unit_t {};
  using height_t = quantity_t<meters_t>;

  inline constexpr depth_t operator"" _gcm(const long double value){return depth_t{value};}
  inline constexpr depth_t operator"" _gcm(const unsigned long long value){return depth_t{static_cast<long double>(value)};}
  inline constexpr height_t operator"" _m(const long double value){return height_t{value};}
  inline constexpr height_t operator"" _cm(const long double value){return height_t{0.01*value};}

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

      m_a = {-186.555305_gcm, -94.919_gcm, 0.61289_gcm, 0.0_gcm, 0.01128292_gcm};
      m_b = {1222.6562_gcm, 1144.9069_gcm, 1305.5948_gcm, 540.1778_gcm, 1_gcm};
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

      if (ilayer < m_nlayers-1) {
        return m_a[ilayer] + m_b[ilayer] * std::exp(-h/m_c[ilayer]);
      } else {
        return m_a[ilayer] - m_b[ilayer]*(h/m_c[ilayer]);
      }
    }

    // * compute height given vertical mass overburden
    height_t
    get_height(
      const depth_t depth
    ) const
    {
      const int ilayer = get_layer_index(depth);

      if (ilayer < m_nlayers-1) {
        return m_c[ilayer] * std::log(m_b[ilayer]/(depth-m_a[ilayer]));
      } else {
        return m_c[ilayer] * ((m_a[ilayer] - depth) / m_b[ilayer]);
      }
    }

    // // * air density as a function of height
    // density_t
    // get_density(
    //   const height_t height
    // ) const
    // {
    //   const int ilayer = get_layer_index(height);

    //   density_t rho = m_b[ilayer]/m_c[ilayer];
    //   if (ilayer < m_nlayers-1) {
    //     rho *= std::exp(-height/m_c[ilayer]);
    //   }

    //   return rho;
    // }

    // // * air density as a function of depth
    // density_t
    // get_density(
    //   const depth_t depth
    // ) const
    // {
    //   const int ilayer = get_layer_index(depth);

    //   if (ilayer < m_nlayers-1) {
    //     return (depth - m_a[ilayer]) / m_c[ilayer];
    //   } else {
    //     return m_b[ilayer]/m_c[ilayer];
    //   }
    // }

  };

}

#endif