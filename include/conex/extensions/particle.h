#ifndef _conex_extensions_particle_h
#define _conex_extensions_particle_h

#include <memory>

#include <conex/extensions/projectile.h>

#include <util/constants.h>
#include <util/vector.h>
#include <util/frame.h>
#include <util/point.h>

namespace conex::extensions {

  class particle;
  using particle_ptr = std::shared_ptr<particle>;

  class particle {

  // * data_t holds the data to read the particle tree
  public: struct data_t {
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
    util::frame_ptr<double> m_frame;
    projectile_ptr m_precursor;

  public:

    // - Constructor
    particle(
      const data_t& tree_data,
      const util::frame_ptr<double>& lab_frame,
      const std::shared_ptr<projectile> precursor
    )
    : m_data(tree_data),
      m_frame(lab_frame),
      m_precursor(precursor)
    {}

    // - Direct access to tree data
    const data_t& data() const
    {return m_data;}

    // - Access to the precursor projectile and other data
    const projectile_ptr& get_precursor() const
    {return m_precursor;}

    const util::frame_ptr<double>& get_frame() const
    {return m_frame;}

    util::point_d get_formation_point() const
    {return get_precursor()->get_position();}

    double get_formation_time_s() const
    {return get_precursor()->get_time_s();}

    // - Formatted access to the particle data

    // * get particle momentum in the lab frame [GeV]
    util::vector_d get_momentum() const
    {return{data().Px, data().Py, data().Pz, m_frame};}

    // * get particle energy [GeV]
    const double& get_energy() const
    {return data().Energy;}

    // * get particle mass [GeV]
    const double& get_mass() const
    {return data().mass;}

    // * get particle velocity [m/s]
    util::vector_d get_velocity() const
    {return get_momentum() * (util::constants::c/get_energy());}

    // * get particle id (conex/nexus code)
    const int& get_id() const
    {return data().id;}

    // * interaction counter
    const int& get_interaction_counter() const
    {return data().interactionCounter;}

  };

} // namespace conex::extensions

#endif // _conex_extensions_particle_h