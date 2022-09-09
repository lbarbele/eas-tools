#ifndef _conex_extensions_projectile_h
#define _conex_extensions_projectile_h

#include <cmath>

#include <util/frame.h>
#include <util/point.h>
#include <util/math.h>
#include <util/vector.h>
#include <util/constants.h>
#include <util/rotation_matrix.h>
#include <util/units.h>

namespace conex::extensions {

  class projectile;
  using projectile_ptr = std::shared_ptr<projectile>;

  class projectile {
  
  // * data_t holds the data to read the projectile tree
  public: struct data_t {
      double Energy = 0;         // dptl(1)
      double Px = 0;             // dptl(2)
      double Py = 0;             // dptl(3)
      double Pz = 0;             // dptl(4)
      double mass = 0;           // dptl(5)
      double x = 0;              // dptl(6)
      double y = 0;              // dptl(7)
      double height = 0;         // dptl(8)
      double time = 0;           // dptl(9)
      double id = 0;             // dptl(10)
      double weight = 0;         // dptl(11)
      double generation = 0;     // dptl(12)
      double slantTraversed = 0; // dptl(13)
      double xShower = 0;        // dptl(14)
      double yShower = 0;        // dptl(15)
      
      int interactionCounter = 0;
      
      double c0s = 0;
      double c0xs = 0;
      double s0s = 0;
      double s0xs = 0;
      
      double slantToImpact = 0;
    };

  private:
    data_t m_data;
    util::frame_ptr<double> m_frame;
    util::frame_ptr<double> m_lab_frame;
    util::point_d m_position;
    util::vector_d m_momentum;

  public:

    // - Constructor
    projectile(
      const data_t& tree_data
    )
    : m_data(tree_data)
    {
      namespace ct = util::constants;
      using util::axis;

      using namespace units::literals;

      // auxiliar quantities
      const double r_obs = std::hypot(data().x, data().y);
      const double dist_center = data().height + ct::earth_radius;
      const double z = std::sqrt((dist_center+r_obs)*(dist_center-r_obs)) - ct::earth_radius;

      // the particle frame
      const auto phi = util::math::atan2(data().y, data().x);
      const auto theta_ea = util::math::asin(r_obs / dist_center);
      m_frame = util::frame<double>::create((axis::x, 1_pi - theta_ea)*(axis::z, phi - 0.5_pi), util::frame<double>::conex_observer);

      // the lab frame (in which secondaries are produced)
      const auto phi_lab = util::math::atan2(data().s0xs, data().c0xs);
      const auto theta_lab = util::math::atan2(data().s0s, data().c0s);
      m_lab_frame = util::frame<double>::create((axis::x, -theta_lab)*(axis::z, -phi_lab), m_frame);

      // position vector
      m_position = util::point_d(data().x, data().y, z, util::frame<double>::conex_observer);

      // momentum vector
      m_momentum = util::vector_d(data().Px, data().Py, data().Pz, m_frame);
    }

    // - Direct access to tree data
    const data_t& data() const
    {return m_data;}

    // - Formatted access to the projectile data

    // * access to the conex "particle" frame
    const util::frame_ptr<double>& get_frame() const
    {return m_frame;}

    // * access to the conex "lab" frame, in which secondary momenta is defined
    const util::frame_ptr<double>& get_lab_frame() const
    {return m_lab_frame;}

    // * get position [m]
    util::point_d get_position() const
    {return m_position;}

    // * get momentum [GeV]
    util::vector_d get_momentum() const
    {return m_momentum;}

    // * get projectile time [s]
    double get_time_s() const
    {return data().time/util::constants::c;}

    // * get projectile energy [GeV]
    const double& get_energy() const
    {return data().Energy;}

    // * get projectile energy [GeV]
    const double& get_mass() const
    {return data().mass;}

    // * get projectile height [m]
    const double& get_height() const
    {return data().height;}

    // * get projectile weight, if thinning is enabled
    const double& get_weight() const
    {return data().weight;}

    // * get projectile id as defined in conex (nexus)
    int get_id() const
    {return data().id;}

    // * projectile generation
    int get_generation() const
    {return data().generation;}

    // * accumulated slant depth traversed by projectile [g/cm^2]
    const double& get_slant_depth() const
    {return data().slantTraversed;}

    // * slant distance to the impact point [m]
    const double& get_distance_to_impact() const
    {return data().slantToImpact;}

    // * projectile x position on the shower plane [m]
    const double& get_xshower() const
    {return data().xShower;}

    // * projectile y position on the shower plane [m]
    const double& get_yshower() const
    {return data().yShower;}

    // * interaction counter
    const int& get_interaction_counter() const
    {return data().interactionCounter;}
  };

}

#endif