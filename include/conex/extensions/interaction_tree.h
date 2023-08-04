#ifndef _conex_extensions_interaction_tree_h
#define _conex_extensions_interaction_tree_h

#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>

#include <conex/extensions/event.h>
#include <conex/extensions/interaction.h>
#include <conex/extensions/particle.h>
#include <conex/extensions/projectile.h>

#include <models/atmosphere.h>

#include <util/vector.h>
#include <util/point.h>

namespace conex::extensions {

  class interaction_tree;
  using interaction_tree_ptr = std::shared_ptr<interaction_tree>;

  class interaction_tree {
  private:
    interaction_ptr m_interaction;
    std::vector<particle_ptr> m_products;
    std::list<interaction_tree_ptr> m_secondary_interactions;
    interaction_tree_ptr m_precursor_interaction;

    units::energy_t m_energy_threshold;
    units::energy_t m_lost_energy;
    int m_generation;

    interaction_tree() {}

    struct transformation_frames {
      util::frame_ptr shower_old = nullptr;
      util::frame_ptr shower_new = nullptr;
      util::frame_ptr observer_new = nullptr;
    };

    interaction_tree_ptr
    do_transform(
      const transformation_frames& frames,
      const util::point_t<units::length_t> initial_position,
      const units::time_t initial_time,
      const units::depth_t traversed_mass
    ) const;

  public:

    // - Transformation

    interaction_tree_ptr
    transform(
      const units::angle_t new_zenith,
      const units::length_t new_observation_level
    ) const;

    // - Static factory methods

    static interaction_tree_ptr create(
      const event& evt,
      const units::energy_t energy_threshold,
      const bool verbose = false
    );

    static interaction_tree_ptr create(
      const interaction_ptr& source,
      std::list<interaction_ptr>& others,
      const units::energy_t energy_threshold,
      const int generation,
      const bool verbose = false
    );

    // - Getters

    const interaction_ptr& get_interaction() const
    {return m_interaction;}

    const projectile_ptr& get_projectile() const
    {return get_interaction()->get_projectile();}

    const std::vector<particle_ptr>& get_final_products() const
    {return m_products;}

    const std::list<interaction_tree_ptr>& get_secondary_interactions() const
    {return m_secondary_interactions;}

    int get_generation() const
    {return m_generation;}

    const interaction_tree_ptr& get_precursor() const
    {return m_precursor_interaction;}

    const units::energy_t& get_energy_loss() const
    {return m_lost_energy;}

    // - Methods

    template <class FcnT>
    void apply_recursive(FcnT fcn) {
      fcn(*this);

      for (interaction_tree_ptr& sub_tree : m_secondary_interactions) {
        sub_tree->apply_recursive(fcn);
      }
    }

  };

}

#endif