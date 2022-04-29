#ifndef _conex_extensions_projectile_h
#define _conex_extensions_projectile_h

#include <cmath>

#include <util/vector.h>

namespace conex::extensions {

  class projectile {
    friend class event;

  public:
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
  
  public:

    util::vector get_position() const;
    util::vector get_momentum() const;

    const double& get_energy() const
    {return Energy;}

    const double& get_px() const
    {return Px;}
    const double& get_py() const
    {return Py;}
    const double& get_pz() const
    {return Pz;}

    const double& get_mass() const
    {return mass;}

    const double& get_x() const
    {return x;}
    const double& get_y() const
    {return y;}
    const double& get_height() const
    {return height;}

    const double& get_time() const
    {return time;}

    int get_id() const
    {return id;}

    const double& get_weight() const
    {return weight;}

    int get_generation() const
    {return generation;}

    const double& get_slant_depth() const
    {return slantTraversed;}

    const int& get_interaction_counter() const
    {return interactionCounter;}

    const double& get_distance_to_impact() const
    {return slantToImpact;}

    const double& get_xshower() const
    {return xShower;}
    const double& get_yshower() const
    {return yShower;}

    // extra info
    double get_r() const
    {return std::hypot(x,y);}

    double get_phip() const
    {return std::atan2(y,x);}

    double get_thetap() const
    {return std::asin(get_r()/(height + 6371315.));}
  };

} // namespace conex::extensions

#endif // _conex_extensions_projectile_h