#include <iostream>
#include <string>
#include <memory>

#include <conex/extensions/file.h>
#include <conex/extensions/event.h>

namespace conex::extensions {

  file::file(
    const std::string& fname
  ) :
    TFile(fname.c_str(), "read"),
    m_event_count(0)
  {
    // check if file is open
    if (!is_open()) {
      std::cerr << "unable to open conex particle file " << fname << std::endl;
      return;
    }

    // get the number of events (showers) in this file
    while(true) {
      auto suffix = std::to_string(m_event_count);

      std::string particle_tree_name = "ParticleList_Event" + suffix;
      std::string projectile_tree_name = "Projectile_Event" + suffix;
      std::string interaction_tree_name = "InteractionList_Event" + suffix;
      std::string seed_tree_name = "Seeds_Event" + suffix;

      std::unique_ptr<TTree> particle_tree(Get<TTree>(particle_tree_name.c_str()));
      std::unique_ptr<TTree> projectile_tree(Get<TTree>(projectile_tree_name.c_str()));
      std::unique_ptr<TTree> interaction_tree(Get<TTree>(interaction_tree_name.c_str()));
      std::unique_ptr<TTree> seed_tree(Get<TTree>(seed_tree_name.c_str()));

      if (particle_tree && projectile_tree && interaction_tree && seed_tree) {
        ++m_event_count;
      } else {
        break;
      }
    }

  }

  event
  file::get_event(
    const size_t pos,
    const double threshold,
    const bool check
  )
  {
    auto suffix = std::to_string(pos);

    std::string particle_tree_name = "ParticleList_Event" + suffix;
    std::string projectile_tree_name = "Projectile_Event" + suffix;
    std::string interaction_tree_name = "InteractionList_Event" + suffix;
    std::string seed_tree_name = "Seeds_Event" + suffix;

    auto particle_tree = Get<TTree>(particle_tree_name.c_str());
    auto projectile_tree = Get<TTree>(projectile_tree_name.c_str());
    auto interaction_tree = Get<TTree>(interaction_tree_name.c_str());
    auto seed_tree = Get<TTree>(seed_tree_name.c_str());

    // the event object is responsible to delete the trees
    return event(particle_tree, projectile_tree, interaction_tree, seed_tree, threshold, check);
  }

} // namespace conex::extensions