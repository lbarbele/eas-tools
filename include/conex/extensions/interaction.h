#ifndef _conex_extensions_interaction_h
#define _conex_extensions_interaction_h

#include <vector>
#include <cstddef>
#include <memory>

#include <conex/extensions/projectile.h>
#include <conex/extensions/particle.h>

#include <util/frame.h>

namespace conex::extensions {

  class interaction {

  // * data_t holds the data to read the interaction and seeds trees
  public: struct data_t {
      int idProj = 0;   // id of projectile (always 1120 if nucleus)
      int idTarg = 0;   // id of target
      int mult = 0;     // multiplicity of secondary particles
      double eProj = 0; // projectile energy in lab frame (energy/nucleon if nucleus)
      double eCMS = 0;  // total energy in CMS
      double eProd = 0; // total energy in lab frame (including target rest mass)
      int interactionCounter = 0;
      std::array<int, 3> seeds;
    };

  private:
    std::shared_ptr<projectile> m_projectile;
    std::vector<particle> m_secondaries;
    data_t m_data;

  public:

    // - Constructor and setters

    // * build from tree data
    interaction(
      const data_t& tree_data
    )
    : m_data(tree_data)
    {}

    // * set interaction projectile
    void set_projectile(const projectile::data_t& proj_data)
    {m_projectile = std::make_shared<projectile>(proj_data);}

    // * add secondary particle to the list
    particle& add_particle(const particle::data_t& part_data)
    {return m_secondaries.emplace_back(part_data, get_lab_frame(), m_projectile);}

    // - Direct access to tree data
    const data_t& data() const
    {return m_data;}

    // - Access to conex frames

    const util::frame_ptr& get_frame() const
    {return get_projectile().get_frame();}

    const util::frame_ptr& get_lab_frame() const
    {return get_projectile().get_lab_frame();}

    // - Formatted access to the tree data

    const int& get_multiplicity() const
    {return data().mult;}
    
    const int& get_projectile_id() const
    {return data().idProj;}

    const int& get_target_id() const
    {return data().idTarg;}

    const double& get_cms_energy() const
    {return data().eCMS;}

    const double& get_lab_energy() const
    {return data().eProd;}

    const double& get_proj_energy() const
    {return data().eProj;}

    const int& get_interaction_counter() const
    {return data().interactionCounter;}

    const std::array<int, 3>& get_seeds() const
    {return data().seeds;}

    // - Access to projectile/secondary particles

    const projectile& get_projectile() const
    {return *m_projectile;}

    const particle& get_secondary(size_t pos) const
    {return m_secondaries[pos];}

    // - Particle iterators

    // * standard iterators
    auto begin() {return m_secondaries.begin();}
    auto end() {return m_secondaries.end();}
    
    // * const iterators
    auto begin() const {return m_secondaries.begin();}
    auto cbegin() const {return m_secondaries.cbegin();}
    auto end() const {return m_secondaries.end();}
    auto cend() const {return m_secondaries.cend();}

    // * reverse iterators
    auto rbegin() {return m_secondaries.rbegin();}
    auto rend() {return m_secondaries.rend();}

    // * const reverse iterators
    auto rbegin() const {return m_secondaries.rbegin();}
    auto crbegin() const {return m_secondaries.crbegin();}
    auto rend() const {return m_secondaries.rend();}
    auto crend() const {return m_secondaries.crend();}
  };

} // namespace conex::extensions

#endif // _conex_extensions_interaction_h