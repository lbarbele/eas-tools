#include <conex/extensions/interaction_tree.h>
#include <util/vector.h>
#include <util/units.h>

// - Helper functions

using namespace units::literals;

namespace {
  bool
  match(
    const conex::extensions::particle_ptr& part,
    const conex::extensions::projectile_ptr& proj,
    const bool verbose = false
  )
  {
    // * check for matching IDs
    // both particle and projectile must have the same particle ID
    if (part->get_id() != proj->get_id()) {
      return false;
    }

    // * check for matching energy
    // requirements are that the projectile energy does not exceed the particle
    // energy (meaning the particle cannot gain energy during propagation) and
    // that relative energy lost is below 0.001%
    const double edev = proj->get_energy()/part->get_energy() - 1;
    if (edev > 1e-10 || edev < -1e-5) {
      if (verbose) std::cerr << "(" << edev << ") edev ";
      return false;
    }

    // * check for matching momentum direction
    // it is required that the angle between particle and projectile momentum
    // is below 10^-7 rad
    const auto angle = util::angle(part->get_momentum(), proj->get_momentum());
    if (angle >= 1e-7_rad) {
      if (verbose) std::cerr << "angle ";
      return false;
    }

    // * check projectile position
    // it is required that, if particle is propagated up to the projectile interaction
    // time, its distance to the interaction point is smaller than 1 micrometer
    const units::time_t deltaTime = proj->get_time() - part->get_formation_time();
    const util::vector_t<units::length_t> posDif =
      part->get_formation_point()  /* initial position */ +
      deltaTime * part->get_velocity() /* displacement */ -
      proj->get_position()    /* actual final position */ ;

    const units::length_t posDifCut = 
      (util::math::abs(proj->get_id()) == 12 || proj->get_id() == 10) ?
      units::length_t(1_m) : units::length_t(1_um);

    if (posDif.norm() > posDifCut) {
      if (verbose) std::cerr << "(" << posDif.norm() << ") " << "position ";
      return false;
    }

    return true;
  }
}

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

    if (verbose) {
      std::cerr
        << "creating interaction tree from event with "
        << evt.get_n_interactions()
        << " interactions\n";
    }

    // get the other parameters
    interaction_ptr source = evt.get_interaction(0);
    const int generation = 0;

    if (verbose) {
      std::cerr
        << "primary interaction data:\n"
        << "+ projectile: " << source->get_projectile()->get_id() << std::endl
        << "+ energy: " << source->get_projectile()->get_energy() << std::endl
        << "+ momentum: " << source->get_projectile()->get_momentum().on_frame(util::frame::conex_observer) << std::endl
        << "+ secondary multiplicity: " << source->get_multiplicity() << std::endl;
    }

    // call overload accepting interaction list
    return create(source, others, energy_threshold, generation, verbose);
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
    // remove the source interaction from the list, if present
    others.remove(source);

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
    tree->m_generation = generation;

    // loop over particles in the current interaction
    for (const particle_ptr& current_particle : *source) {

      // dismiss particles with low kinectic energy (say, 10 MeV)
      if (current_particle->get_energy() - (current_particle->get_mass()*1_c*1_c) < 10_MeV) {
        continue;
      }

      // do not check particles with energy below threshold, simply add them to the stack.
      // also, if there are no other interactions to analyze, simply add the particle
      // to the stack
      if (current_particle->get_energy() < energy_threshold || others.empty()) {
        tree->m_products.emplace_back(current_particle);
        continue;
      }

      // pointer to the matching interaction, if any
      interaction_ptr matching_interaction = nullptr;

      if (verbose) {
        std::cerr
          << std::endl
          << std::left
          << "+ searching "
          << "(id) " << std::setw(5) << current_particle->get_id()
          << " (E) " << std::setw(15) << current_particle->get_energy()
          << " (p) " << std::setw(15) << current_particle->get_momentum().on_frame(util::frame::conex_observer)
          << std::right << std::endl;
      }

      // iterate over remaining interactions, searching for a match with the current particle
      for (const interaction_ptr& interaction : others) {

        if (verbose && interaction->get_projectile_id() == current_particle->get_id()) {
          std::cerr
            << std::left
            << ". candidate "
            << "(id) " << std::setw(5) << interaction->get_projectile_id()
            << " (E) " << std::setw(15) << interaction->get_proj_energy()
            << " (p) " << std::setw(15) <<interaction->get_projectile()->get_momentum().on_frame(util::frame::conex_observer)
            << " ... " << std::right;
        }

        // check if this is a match. if so, get a pointer to such interaction in case this is
        // the first match. if it is a match, but not the first one, disambiguate by selecting
        // the interaction that happened first
        if (match(current_particle, interaction->get_projectile(), verbose)) {
          if (!matching_interaction || interaction->get_time() < matching_interaction->get_time()) {
            matching_interaction = interaction;
          }
          if (verbose) {
            std::cerr << "match!" << std::endl;
          }
        } else if (verbose && interaction->get_projectile_id() == current_particle->get_id()) {
          std::cerr << "fail" << std::endl;
        }
      }

      // if there were no matches, add particle to the list of final products, then
      // continue with the next particle
      if (!matching_interaction) {
        // notify if high energy particle did not match any interaction
        std::cerr
          << "! unmatched particle ("
          << "ID: " << current_particle->get_id()
          << ", E: " << current_particle->get_energy()
          << ", E/Ethr: " << current_particle->get_energy()/energy_threshold
          << ")\n";
        tree->m_products.push_back(current_particle);
        continue;
      }

      // add the matching interaction as a sub-interaction_tree relative to this interaction tree
      interaction_tree_ptr subtree = create(
        matching_interaction,
        others,
        energy_threshold,
        generation+1,
        verbose
      );

      // add this interaction as precursor of the subtree
      subtree->m_precursor_interaction = tree;

      // add the subtree as a secondary interaction to this
      tree->m_secondary_interactions.emplace_back(subtree);
    } // loop over particles in the interaction

    // return the newly created tree
    return tree;
  }

}