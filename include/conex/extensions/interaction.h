#ifndef _conex_extensions_interaction_h
#define _conex_extensions_interaction_h

#include <vector>
#include <cstddef>

#include <conex/extensions/projectile.h>
#include <conex/extensions/particle.h>

namespace conex::extensions {

  class interaction {
    friend class event;

  private:
    projectile m_projectile;
    std::vector<particle> m_secondaries;

    // interaction tree
    int idProj = 0;
    int idTarg = 0;
    int mult = 0;
    double eProj = 0;
    double eCMS = 0;
    double eProd = 0; // ? what is this?
    int interactionCounter = 0;

    // seed tree
    int seed1 = 0;
    int seed2 = 0;
    int seed3 = 0;

  public:
    const projectile& get_projectile() const
    {return m_projectile;}

    const particle& get_secondary(size_t pos) const
    {return m_secondaries[pos];}

    size_t get_multiplicity() const
    {return mult;}
    
    int get_projectile_id() const
    {return idProj;}

    int get_target_id() const
    {return idTarg;}

    int get_interaction_counter() const
    {return interactionCounter;}

    double get_cms_energy() const
    {return eCMS;}

    double get_lab_energy() const
    {return eProj;}

    // iteration through secondary particles
    auto begin() const {return m_secondaries.cbegin();}
    auto cbegin() const {return m_secondaries.cbegin();}
    auto end() const {return m_secondaries.cend();}
    auto cend() const {return m_secondaries.cend();}
  };

} // namespace conex::extensions

#endif // _conex_extensions_interaction_h