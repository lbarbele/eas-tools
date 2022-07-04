#ifndef _conex_extensions_interaction_tree_h
#define _conex_extensions_interaction_tree_h

#include <iostream>
#include <list>
#include <vector>
#include <memory>

#include <conex/extensions/event.h>
#include <conex/extensions/interaction.h>
#include <conex/extensions/particle.h>
#include <conex/extensions/projectile.h>

namespace conex::extensions {

  class interaction_tree;
  using interaction_tree_ptr = std::shared_ptr<interaction_tree>;

  class interaction_tree {
  private:
    interaction_ptr m_interaction;
    std::vector<particle_ptr> m_products;
    std::list<interaction_tree_ptr> m_secondary_interactions;
    interaction_tree_ptr m_precursor_interaction;

    double m_energy_threshold;
    int m_generation;

    interaction_tree() {}

  public:

    // - Static creators

    static interaction_tree_ptr create(
      const event& evt,
      const double energy_threshold
    );

    static interaction_tree_ptr create(
      const interaction_ptr& source,
      std::list<interaction_ptr>& others,
      const double energy_threshold,
      const int generation
    );

    // - Getters

    const interaction_ptr& get_interaction() const
    {return m_interaction;}

    const std::vector<particle_ptr>& get_final_products() const
    {return m_products;}

    const std::list<interaction_tree_ptr>& get_secondary_interactions() const
    {return m_secondary_interactions;}

    int get_generation() const
    {return m_generation;}

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