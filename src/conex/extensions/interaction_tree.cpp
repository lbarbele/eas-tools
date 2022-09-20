#include <utility>

#include <conex/extensions/interaction_tree.h>
#include <util/vector.h>
#include <util/units.h>

// - Helper functions

using namespace units::literals;

namespace conex::extensions {

  interaction_tree_ptr
  interaction_tree::create(
    const event& evt,
    const units::energy_t energy_threshold,
    const bool verbose
  )
  {
    // copy the vector of interactions for the event into a list
    std::list<interaction_ptr> others(evt.begin(), evt.end());

    // get the source interaction
    interaction_ptr source = evt.get_interaction(0);

    // remove the source interaction from the list of interactions to analyze
    others.remove(source);

    // start the generation counter
    const int generation = 0;

    // call overload accepting interaction list
    auto tree = create(source, others, energy_threshold, generation, verbose);

    // if verbose, print how many interactions were left unmatched
    if (verbose) {
      std::cerr << "\ndone with " << others.size() << " interactions left\n";
    }

    return tree;
  }

  interaction_tree_ptr
  interaction_tree::create(
    const interaction_ptr& source,
    std::list<interaction_ptr>& others,
    const units::energy_t energy_threshold,
    const int generation,
    const bool verbose
  )
  {
    // printout
    if (verbose) {
      std::cerr << std::endl
        << "creating interaction tree from source:\n"
        << "+ projectile: " << source->get_projectile()->get_id() << std::endl
        << "+ energy: " << source->get_projectile()->get_energy() << std::endl
        << "+ momentum: " << source->get_projectile()->get_momentum().on_frame(util::frame::conex_observer) << std::endl
        << "+ secondary multiplicity: " << source->get_multiplicity() << std::endl;
    }

    // create an empty interaction tree and set its basic properties
    interaction_tree_ptr tree(new interaction_tree);
    tree->m_interaction = source;
    tree->m_energy_threshold = energy_threshold;
    tree->m_lost_energy = units::energy_t(0);
    tree->m_generation = generation;
 
    // create a (empty) list of matching interactions
    std::list<std::pair<particle_ptr, interaction_ptr>> matches;

    // loop over particles in the current interaction
    for (const particle_ptr& current_particle : *source) {

      // dismiss particles with low kinectic energy (say, 10 MeV)
      if (current_particle->get_energy() - (current_particle->get_mass()*1_c*1_c) < 10_MeV) {
        tree->m_lost_energy += current_particle->get_energy();
        continue;
      }

      // do not check particles with energy below threshold, simply add them to the stack.
      // also, if there are no other interactions to analyze, simply add the particle
      // to the stack
      if (current_particle->get_energy() < energy_threshold || others.empty()) {
        tree->m_products.emplace_back(current_particle);
        continue;
      }

      // print the particle whose interaction we are looking for
      if (verbose) {
        std::cerr
          << std::left
          << "- searching "
          << "(id) " << std::setw(5) << current_particle->get_id()
          << " (E) " << std::setw(15) << current_particle->get_energy()
          << " (p) " << std::setw(15) << current_particle->get_momentum().on_frame(util::frame::conex_observer)
          << " among " << others.size() << " interactions ... "
          << std::right;
      }

      // pointer to the matching interaction, if any
      interaction_ptr matching_interaction = nullptr;

      // iterate over remaining interactions, searching for a match with the current particle
      for (const interaction_ptr& interaction : others) {
        if (current_particle->get_unique_id() == interaction->get_projectile()->get_unique_id()) {
          matching_interaction = interaction;
          tree->m_lost_energy += current_particle->get_energy() - interaction->get_projectile()->get_energy();
          break;
        }
      }

      // tell if search succeeded
      if (verbose) {
        std::cerr << (matching_interaction? "ok" : "fail") << std::endl;
      }

      if (!matching_interaction) {
        // no matching interaction

        // notify if high energy particle did not match any interaction because
        // particles with energy above threshold are expected to match
        std::cerr
          << "! unmatched particle ("
          << "ID: " << current_particle->get_id()
          << ", E: " << current_particle->get_energy()
          << ", E/Ethr: " << current_particle->get_energy()/energy_threshold
          << ")\n";

        // add particle to the stack of final products of this interaction
        tree->m_products.push_back(current_particle);

        // go to the next particle
        continue;
      } else {
        // add interaction to the list of matches
        matches.emplace_back(current_particle, matching_interaction);
        // remove interaction from the list to be analyzed
        others.remove(matching_interaction);
      }

    } // loop over particles in the interaction

    // print statistics
    if (verbose) {
      std::cerr << ". " << matches.size() << " particle matches" << std::endl;
      std::cerr << ". " << tree->m_products.size() << " final products" << std::endl;
      std::cerr << ". " << source->get_multiplicity()-matches.size()-tree->m_products.size() << " discarded" << std::endl;
      std::cerr << ". " << others.size() << " interactions remaining" << std::endl;
      std::cerr << ". energy loss was " << tree->m_lost_energy << std::endl;
    }

    // loop over matches and create the subtrees
    for (const auto& [particle, interaction] : matches) {
      // add the matching interaction as a sub-interaction_tree relative to this interaction tree
      interaction_tree_ptr subtree =
        create(interaction, others, energy_threshold, generation+1, verbose);

      // add this interaction as precursor of the subtree
      subtree->m_precursor_interaction = tree;

      // add the subtree as a secondary interaction to this
      tree->m_secondary_interactions.emplace_back(subtree);
    }

    // return the newly created tree
    return tree;
  }

}