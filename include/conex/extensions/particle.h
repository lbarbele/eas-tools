#ifndef _conex_extensions_particle_h
#define _conex_extensions_particle_h

#include <memory>

#include <conex/extensions/projectile.h>

#include <util/constants.h>
#include <util/vector.h>
#include <util/frame.h>
#include <util/point.h>
#include <util/units.h>

namespace conex::extensions {

  // - forward declarations and aliases

  // * the particle class
  class particle;

  // * smart pointer to the particle class
  using particle_ptr = std::shared_ptr<particle>;

  // - implementation of the particle class

  class particle {
  public:
  
    // * data_t holds the data to read the particle tree
    struct data_t {
      double Px = 0;            // xsptl(1,i) ...... x-component of particle momentum 
      double Py = 0;            // xsptl(2,i) ...... y-component of particle momentum 
      double Pz = 0;            // xsptl(3,i) ...... z-component of particle momentum 
      double Energy = 0;        // xsptl(4,i) ...... particle energy 
      double mass = 0;          // xsptl(5,i) ...... particle mass 

      double x = 0;             // xsorptl(1,i) .... x-component of formation point
      double y = 0;             // xsorptl(2,i) .... y-component of formation point
      double z = 0;             // xsorptl(3,i) .... z-component of formation point
      double time = 0;          // xsorptl(4,i) .... formation time

      int id = 0;               // idptlxs(i) ...... particle id

      double t_formation = 0;   // xstivptl(1,i) ... formation time (always in the pp-cms!)
      double t_destruction = 0; // xstivptl(2,i) ... destruction time (always in the pp-cms!)

      int id_origin = 0;        // ityptlxs(i)  .... type of particles origin:

      int id_father = 0;        // iorptlxs(i) ..... particle number of father (if .le. 0 : no father)

      int id_mother = 0;        // jorptlxs(i) ..... particle number of mother (if .le. 0 : no mother)

      int status = 0;           // istptlxs(i) ..... status

      int interactionCounter = 0;
    };
  
  private:
    data_t m_data;
    util::frame_ptr m_frame;
    projectile_ptr m_precursor;

  public:

    // * constructor
    particle(
      const data_t& tree_data,
      const util::frame_ptr& lab_frame,
      const std::shared_ptr<projectile> precursor
    )
    : m_data(tree_data),
      m_frame(lab_frame),
      m_precursor(precursor)
    {}

    // * direct access to raw tree data
    const data_t& data() const
    {return m_data;}

    // * access to the precursor projectile and other data

    const projectile_ptr& get_precursor() const
    {return m_precursor;}

    const util::frame_ptr& get_frame() const
    {return m_frame;}

    util::point_t<units::length_t> get_formation_point() const
    {return get_precursor()->get_position();}

    units::time_t get_formation_time() const
    {return get_precursor()->get_time();}

    // * formatted access to the particle data

    // get particle momentum in the lab frame
    util::vector_t<units::momentum_t> get_momentum() const
    {
      return {
        units::gev_per_c_t<double>(data().Px),
        units::gev_per_c_t<double>(data().Py),
        units::gev_per_c_t<double>(data().Pz),
        m_frame};
    }

    // get particle energy
    units::energy_t get_energy() const
    {return units::gigaelectron_volt_t<double>(data().Energy);}

    // get particle mass
    units::mass_t get_mass() const
    {return units::gev_per_c_squared_t<double>(data().mass);}

    // get particle velocity
    util::vector_t<units::speed_t> get_velocity() const
    {
      using namespace units::literals;
      return get_momentum() * (1_c*1_c/get_energy());
    }

    // get particle id (conex/nexus code)
    const int& get_id() const
    {return data().id;}

    // interaction counter
    const int& get_interaction_counter() const
    {return data().interactionCounter;}

  };

} // namespace conex::extensions

#endif // _conex_extensions_particle_h