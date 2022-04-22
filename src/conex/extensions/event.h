#ifndef _conex_extensions_event_h
#define _conex_extensions_event_h

#include <memory>
#include <vector>

#include <TTree.h>

#include <conex/extensions/interaction.h>
#include <conex/extensions/particle.h>
#include <conex/extensions/projectile.h>

namespace conex::extensions {

  class event {
  public:

    std::unique_ptr<TTree> m_particle_tree;
    std::unique_ptr<TTree> m_projectile_tree;
    std::unique_ptr<TTree> m_interaction_tree;
    std::unique_ptr<TTree> m_seed_tree;

    std::vector<size_t> m_particle_tree_indices;

    interaction m_interaction;
    particle m_particle;

  public:
    event(TTree* particle, TTree* projectile, TTree* interaction, TTree* seed);

    long long get_n_interactions() const;
    interaction get_interaction(size_t pos);
  };

} // namespace conex::extensions

#endif // _conex_extensions_event_h