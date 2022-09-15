#ifndef _conex_extensions_projectile_h
#define _conex_extensions_projectile_h

#include <iostream>
#include <cmath>
#include <memory>

#include <util/constants.h>
#include <util/frame.h>
#include <util/math.h>
#include <util/point.h>
#include <util/rotation_matrix.h>
#include <util/units.h>
#include <util/vector.h>

namespace conex::extensions {

  // - forward declarations and aliases

  // * the projectile class
  class projectile;

  // * smart ptr to the projectile class
  using projectile_ptr = std::shared_ptr<projectile>;

  // - implementation of the projectile class

  class projectile {
  public:
  
    // * data_t holds the data to read the projectile tree
    struct data_t {
      units::gigaelectron_volt_t<double> Energy;                    // dptl(1)
      units::gev_per_c_t<double> Px;                                // dptl(2)
      units::gev_per_c_t<double> Py;                                // dptl(3)
      units::gev_per_c_t<double> Pz;                                // dptl(4)
      units::gev_per_c_squared_t<double> mass;                      // dptl(5)
      units::meter_t<double> x;                                     // dptl(6)
      units::meter_t<double> y;                                     // dptl(7)
      units::meter_t<double> height;                                // dptl(8)
      units::meter_t<double> time;                                  // dptl(9)
      double id = 0;                                                // dptl(10)
      double weight = 0;                                            // dptl(11)
      double generation = 0;                                        // dptl(12)
      units::grams_per_squared_centimeter_t<double> slantTraversed; // dptl(13)
      units::meter_t<double> xShower;                               // dptl(14)
      units::meter_t<double> yShower;                               // dptl(15)
      
      int interactionCounter = 0;
      
      double c0s = 0;
      double c0xs = 0;
      double s0s = 0;
      double s0xs = 0;
      
      double slantToImpact = 0;
    };

  private:
    data_t m_data;
    util::frame_ptr m_frame;
    util::frame_ptr m_lab_frame;
    util::point_t<units::length_t> m_position;
    util::vector_t<units::momentum_t> m_momentum;

  public:

    // * constructor from raw tree data

    projectile(const data_t& tree_data) : m_data(tree_data)
    {
      using namespace units::literals;
 
      const units::length_t dist_center = data().height + util::constants::earth_radius;
      const units::length_t r = util::math::hypot(data().x, data().y);
      const units::length_t z = util::math::sqrt((dist_center + r) * (dist_center - r)) - util::constants::earth_radius;

      // the position vector
      m_position = util::point_t<units::length_t>(data().x, data().y, z, util::frame::conex_observer);

      // the particle frame
      const units::angle_t phi_pos = util::math::atan2(data().y, data().x);
      const units::angle_t theta_ea = util::math::asin(r / dist_center);
      const auto particle_frame_rot = (util::axis::x, 1_pi - theta_ea) * (util::axis::z, phi_pos - 0.5_pi);
      m_frame = util::frame::create(particle_frame_rot, m_position, util::frame::conex_observer);

      // the lab frame (in which secondaries are produced)
      const units::angle_t phi_lab = util::math::atan2(data().s0xs, data().c0xs);
      const units::angle_t tht_lab = util::math::atan2(data().s0s, data().c0s);
      const auto lab_frame_rot = (util::axis::x, -tht_lab)*(util::axis::z, -phi_lab);
      m_lab_frame = util::frame::create(lab_frame_rot, m_frame);

      // momentum vector
      m_momentum.set_frame(m_frame);
      m_momentum[0] = data().Px;
      m_momentum[1] = data().Py;
      m_momentum[2] = data().Pz;
    }

    // * direct access to raw tree data

    const data_t& data() const
    {return m_data;}

    // * formatted access to the projectile data

    // access to the conex "particle" frame
    const util::frame_ptr& get_frame() const
    {return m_frame;}

    // access to the conex "lab" frame, in which secondary momenta is defined
    const util::frame_ptr& get_lab_frame() const
    {return m_lab_frame;}

    // get position
    util::point_t<units::length_t> get_position() const
    {return m_position;}

    // get momentum
    util::vector_t<units::momentum_t> get_momentum() const
    {return m_momentum;}

    // get projectile time
    units::time_t get_time() const
    {return data().time/util::constants::speed_of_light;}

    // get projectile energy
    units::energy_t get_energy() const
    {return data().Energy;}

    // get projectile mass [GeV / c^2]
    units::mass_t get_mass() const
    {return data().mass;}

    // get projectile velocity
    util::vector_t<units::speed_t> get_velocity() const
    {
      using namespace units::literals;

      const util::vector_d beta = get_id() == 10?
        get_momentum().get_normalized(1) :
        get_momentum() * 1_c / get_energy();

      return 1_c * beta;
    }

    // get projectile height [m]
    units::length_t get_height() const
    {return units::meter_t<double>(data().height);}

    // get projectile weight, if thinning is enabled
    const double& get_weight() const
    {return data().weight;}

    // get projectile id as defined in conex (nexus)
    int get_id() const
    {return data().id;}

    // projectile generation
    int get_generation() const
    {return data().generation;}

    // accumulated slant depth traversed by projectile [g/cm^2]
    units::depth_t get_slant_depth() const
    {return data().slantTraversed;}

    // slant distance to the impact point [m]
    units::length_t get_distance_to_impact() const
    {return units::meter_t<double>(data().slantToImpact);}

    // projectile x position on the shower plane [m]
    units::length_t get_xshower() const
    {return data().xShower;}

    // projectile y position on the shower plane [m]
    units::length_t get_yshower() const
    {return data().yShower;}

    // interaction counter
    const int& get_interaction_counter() const
    {return data().interactionCounter;}
  };

}

#endif