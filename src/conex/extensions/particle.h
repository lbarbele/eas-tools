#ifndef _conex_extensions_particle_h
#define _conex_extensions_particle_h

#include <util/vector.h>

namespace conex::extensions {

  class particle {
    friend class event;

  private:
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

  public:
    util::vector get_momentum() const;

    double get_energy() const
    {return Energy;}

    double get_mass() const
    {return mass;}

    double get_formation_time() const
    {return t_formation;}

    double get_destruction_time() const
    {return t_destruction;}

    int get_id() const
    {return id;}

    int get_interaction_counter() const
    {return interactionCounter;}
  };

} // namespace conex::extensions

#endif // _conex_extensions_particle_h