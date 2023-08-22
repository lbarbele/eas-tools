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

  interaction_tree_ptr
  interaction_tree::transform(
    const units::angle_t new_zenith,
    const units::length_t new_observation_level
  ) const
  {
    // * the transformation is always applied to the head of the tree

    if (get_precursor() != nullptr) {
      return get_precursor()->transform(new_zenith, new_observation_level);
    }

    // * helpers

    using namespace units::literals;
    const auto rea = util::constants::earth_radius;
    const auto atm = models::atmosphere::us_standard();
    const auto conex_observer = util::frame::conex_observer;

    // * a container for the frames

    transformation_frames frames;

    // * define the old shower frame

    const auto old_projectile = get_interaction()->get_projectile();
    const auto old_momentum = old_projectile->get_momentum();
    const auto old_axis = -old_momentum.get_normalized(1.).on_frame(conex_observer);

    const auto old_theta = util::math::acos(old_axis.z());
    const auto old_phi = util::math::atan2(old_axis.y(), old_axis.x());

    // the shower frame is defined such that the direction of the primary-particle momentum
    // coincides with the z axis (points downwards). the y axis is such that its projection
    // on the ground coincides with the projection of the primary particle momentum.
    // the origin of this frame is also displaced to an altitude given by the observation
    // level.
    frames.shower_old = util::frame::create(
      (util::axis::x, 1_pi - old_theta)* // 2: align z axis with the shower axis such that z points downwards
      (util::axis::z, old_phi - 0.5_pi), // 1: align y axis with projection of shower axis on ground
      conex_observer
    );

    // * define the new observer frame

    const auto new_obs_point = util::point_t<units::length_t>(0_m, 0_m, new_observation_level, conex_observer);
    frames.observer_new = util::frame::create(new_obs_point, conex_observer);

    // * define the new shower frame

    // new angles of the shower axis and sines/cosines
    const auto new_theta = new_zenith;
    const auto new_phi = old_phi;

    const auto new_sin_theta = util::math::sin(new_theta);
    const auto new_cos_theta = util::math::cos(new_theta);

    frames.shower_new = util::frame::create(
      (util::axis::x, 1_pi - new_theta)* // 2: align z axis with the shower axis such that z points downwards
      (util::axis::z, new_phi - 0.5_pi), // 1: align y axis with projection of shower axis on ground
      frames.observer_new
    );

    // * define the new initial position

    // distance from earth center to the border of the atmosphere
    const auto rmax = rea + atm.max_height();
    // distance from earth center to the observation level
    const auto robs = rea + new_observation_level;
    // helper variables
    const auto xobs = robs * new_cos_theta;
    const auto yobs = robs * new_sin_theta;
    // distance from observation point to atm. border along new shower axis
    const auto new_max_distance = util::math::sqrt((rmax + yobs)*(rmax - yobs)) - xobs;
    // the new initial position is the intercept between shower axis and atmosphere boundary
    const auto new_initial_position = util::point_t<units::length_t>(0_m, 0_m, -new_max_distance, frames.shower_new);

    // * define the new time frame

    const auto new_initial_time = -new_max_distance/1_c;

    // * apply transformation

    const auto traversed_mass = old_projectile->get_slant_depth();

    return do_transform(frames, new_initial_position, new_initial_time, traversed_mass);
  }

  interaction_tree_ptr
  interaction_tree::do_transform(
    const transformation_frames& frames,
    const util::point_t<units::length_t> initial_position,
    const units::time_t initial_time,
    const units::depth_t traversed_mass
  ) const
  {
    // * helpers

    using namespace units::literals;
    const auto atm = models::atmosphere::us_standard();
    const auto conex_observer = util::frame::conex_observer;
    const auto earth_radius = util::constants::earth_radius;
    const auto earth_center = util::point_t<units::length_t>(0_m, 0_m, -earth_radius, conex_observer);
    const auto observation_level = (
      util::point_t<units::length_t>(0_m, 0_m, 0_m, frames.observer_new) -
      util::point_t<units::length_t>(0_m, 0_m, 0_m, util::frame::conex_observer)
    ).norm();

    // * redefine projectile momentum

    // read original momentum in the old shower frame
    const auto old_projectile = get_projectile();
    const auto old_momentum = old_projectile->get_momentum().on_frame(frames.shower_old);

    // define new momentum on the new shower frame
    auto new_momentum = util::vector_t<units::momentum_t>(
      old_momentum.x(),
      old_momentum.y(),
      old_momentum.z(),
      frames.shower_new
    );

    // * move particle from initial position

    const auto final_position = atm
      .transport(initial_position, new_momentum, traversed_mass, observation_level)
      .on_frame(util::frame::conex_observer);

    // new particle height (above sea level!)
    const auto new_height = atm.get_height(final_position);

    // flag indicating the projectile has hit the ground
    const bool particle_hits_ground = new_height - observation_level < 1_nm;

    // new particle time
    const auto displacement = (final_position - initial_position).norm();
    const auto velocity = old_projectile->get_velocity().norm();
    const auto final_time = initial_time + displacement/velocity;

    // * compute particle frame at the final position

    // spherical coordinates of final position with origin at the earth center
    const auto rad = new_height + earth_radius;
    const auto rxy = util::math::hypot(final_position.y(), final_position.x());
    const auto phi = util::math::atan2(final_position.y(), final_position.x());
    const auto tea = util::math::asin(rxy/rad);

    // define the new particle frame, as in CONEX
    const auto particle_frame_rot = (util::axis::x, 1_pi-tea)*(util::axis::z, phi-0.5_pi);
    const auto new_particle_frame = util::frame::create(particle_frame_rot, final_position, util::frame::conex_observer);

    const auto ax1 = util::vector_d(0, 0, 1, get_projectile()->get_frame());
    const auto ax2 = util::vector_d(0, 0, 1, new_particle_frame);

    // * compute new slant distance to impact and new slant depth

    const auto new_shower_axis = util::vector_d(0, 0, -1, frames.shower_new);
    const auto new_slant_distance = -final_position.on_frame(frames.shower_new).z();
    const auto new_slant_depth = atm.get_slant_depth(new_shower_axis, new_slant_distance, observation_level);

    // * compute new euler angles (sines and cosines) for the lab frame

    const auto p = new_momentum.on_frame(new_particle_frame);

    const auto new_phi_lab = util::math::atan2(p.x(), p.y()); // x, y swapped (see sbrt CXDEFROT)
    const auto new_c0xs = util::math::cos(new_phi_lab);
    const auto new_s0xs = util::math::sin(new_phi_lab);
    const auto new_c0s = p.z()/p.norm();
    const auto new_s0s = util::math::sqrt((1 + new_c0s)*(1 - new_c0s));

    // * redefine the projectile data structure

    const auto new_projectile_data = projectile::data_t{
      // the particle momentum defined on the particle frame
      .Px = new_momentum.on_frame(new_particle_frame).x(),
      .Py = new_momentum.on_frame(new_particle_frame).y(),
      .Pz = new_momentum.on_frame(new_particle_frame).z(),
      // energy and mass are unchanged
      .Energy = old_projectile->get_energy(),
      .mass   = old_projectile->get_mass(),
      // the horizontal cartesian coordinates on the conex observer frame
      .x = final_position.on_frame(util::frame::conex_observer).x(),
      .y = final_position.on_frame(util::frame::conex_observer).y(),
      // height is always measured above sea level, independent of the observation level
      .height = new_height,
      // accumulated time in meters (actually t*c)
      .time = final_time * 1_c,
      // weight, id, and generation are unchanged
      .id         = double(old_projectile->get_id()),
      .weight     = old_projectile->get_weight(),
      .generation = double(old_projectile->get_generation()),
      // slant depth computed along shower axis
      .slantTraversed = new_slant_depth,
      // x, y coordinates in the new shower frame
      .xShower = final_position.on_frame(frames.shower_new).x(),
      .yShower = final_position.on_frame(frames.shower_new).y(),
      // slant distance to impact is -z coordinate in the shower frame
      .slantToImpact = new_slant_distance,
      // unique id and interaction couter are unchanged
      .uniqueParticleId   = double(old_projectile->get_unique_id()),
      .interactionCounter = old_projectile->get_interaction_counter(),
      // rotations specifying the interaction (lab) frame
      .c0s  = new_c0s,
      .c0xs = new_c0xs,
      .s0s  = new_s0s,
      .s0xs = new_s0xs,
    };

    // * create a new interaction object, still without secondary particles

    const auto new_interaction_data = interaction::data_t{
      .idProj             = get_interaction()->data().idProj,
      .idTarg             = get_interaction()->data().idTarg,
      .mult               = particle_hits_ground ? 1 : get_interaction()->data().mult,
      .eProj              = get_interaction()->data().eProj,
      .eCMS               = get_interaction()->data().eCMS,
      .eProd              = get_interaction()->data().eProd,
      .interactionCounter = get_interaction()->data().interactionCounter,
      .seeds              = get_interaction()->data().seeds
    };

    const auto new_interaction = std::make_shared<interaction>(new_interaction_data);
    
    // set the projectile
    new_interaction->set_projectile(new_projectile_data);

    // * create a new tree node

    // create empty tree
    const auto new_tree = interaction_tree_ptr(new interaction_tree);

    // copy tree properties
    new_tree->m_energy_threshold = m_energy_threshold;
    new_tree->m_lost_energy = m_lost_energy;
    new_tree->m_generation = m_generation;

    // set interaction data
    new_tree->m_interaction = new_interaction;

    if (particle_hits_ground) {
    // * if projectile hits the observation level

      // add the projectile as the only secondary particle emerging from the interaction
      const auto new_sec_particle = new_interaction->add_particle({
        .Px = units::gev_per_c_t<double>(0.),
        .Py = units::gev_per_c_t<double>(0.),
        .Pz = new_interaction->get_projectile()->get_momentum().norm(),
        .Energy = new_projectile_data.Energy,
        .mass = new_projectile_data.mass,
        .x = 0_m,
        .y = 0_m,
        .z = 0_m,
        .time = 0_s,
        .id = new_projectile_data.id,
        .t_formation = 0_s,
        .t_destruction = 0_s,
        .id_origin = 0,
        .id_father = 0,
        .id_mother = 0,
        .status = 0,
        .uniqueParticleId = new_projectile_data.uniqueParticleId,
        .interactionCounter = new_projectile_data.interactionCounter
      });

      // the projectile is also a final product on the tree node
      new_tree->m_products.push_back(new_sec_particle);

      // note: no subtrees are added!

    } else {
    // * if projectile interacts

      // add all secondary particles to the new interaction
      for (const auto& secondary_particle : *get_interaction()) {
        new_interaction->add_particle(secondary_particle->data());
      }

      // set final products (find by matching indices) of this interaction tree node
      for (const auto& ptc : *new_interaction) {
        bool is_final = false;

        for (const auto& prod : get_final_products()) {
          if (prod->get_unique_id() == ptc->get_unique_id()) {
            is_final = true;
            break;
          }
        }

        if (is_final) {
           new_tree->m_products.push_back(ptc);
        }
      }

      // apply transformation to secondary interactions and connect them to the new tree
      for (const auto& old_subtree : get_secondary_interactions()) {

        // get mass traversed between this interaction and the next (as computed in CONEX)
        const auto dz = atm.get_traversed_mass(
          old_projectile->get_position(),
          old_subtree->get_projectile()->get_position()
        );

        // create a new subtree by transforming the old subtree
        const auto new_subtree = old_subtree->do_transform(frames, final_position, final_time, dz);

        // connect the new subtree to the new tree
        new_subtree->m_precursor_interaction = new_tree;
        new_tree->m_secondary_interactions.push_back(new_subtree);
      }
    }

    return new_tree;
  }

}