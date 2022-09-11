#include <stdexcept>
#include <sstream>

#include <conex/extensions/event.h>
#include <conex/extensions/interaction.h>
#include <conex/extensions/projectile.h>
#include <conex/extensions/particle.h>

#include <util/units.h>
#include <util/vector.h>

#include <iostream>
#include <iomanip>

namespace conex::extensions {

  event::event(
    TTree* particle_tree,
    TTree* projectile_tree,
    TTree* interaction_tree,
    TTree* seed_tree,
    const double threshold,
    const bool check
  )
  {
    using namespace units::literals;

    // * check trees
    if (!particle_tree || !projectile_tree || !interaction_tree || !seed_tree) {
      throw("conex::extensions::event::event(): bad tree pointers");
    }

    // * interaction/seeds tree
    interaction::data_t interaction_data;

    interaction_tree->SetBranchAddress("idProj",&interaction_data.idProj);
    interaction_tree->SetBranchAddress("idTarg",&interaction_data.idTarg);
    interaction_tree->SetBranchAddress("mult",&interaction_data.mult);
    interaction_tree->SetBranchAddress("eProj",&interaction_data.eProj);
    interaction_tree->SetBranchAddress("eCMS",&interaction_data.eCMS);
    interaction_tree->SetBranchAddress("eProd",&interaction_data.eProd);
    interaction_tree->SetBranchAddress("interactionCounter",&interaction_data.interactionCounter);

    seed_tree->SetBranchAddress("seed1a",&interaction_data.seeds[0]);
    seed_tree->SetBranchAddress("seed2a",&interaction_data.seeds[1]);
    seed_tree->SetBranchAddress("seed3a",&interaction_data.seeds[2]);

    // * projectile tree
    projectile::data_t projectile_data;

    projectile_tree->SetBranchAddress("Energy",&projectile_data.Energy);
    projectile_tree->SetBranchAddress("Px",&projectile_data.Px);
    projectile_tree->SetBranchAddress("Py",&projectile_data.Py);
    projectile_tree->SetBranchAddress("Pz",&projectile_data.Pz);
    projectile_tree->SetBranchAddress("mass",&projectile_data.mass);
    projectile_tree->SetBranchAddress("x",&projectile_data.x);
    projectile_tree->SetBranchAddress("y",&projectile_data.y);
    projectile_tree->SetBranchAddress("time",&projectile_data.time);
    projectile_tree->SetBranchAddress("id",&projectile_data.id);
    projectile_tree->SetBranchAddress("interactionCounter",&projectile_data.interactionCounter);
    projectile_tree->SetBranchAddress("c0s",&projectile_data.c0s);
    projectile_tree->SetBranchAddress("c0xs",&projectile_data.c0xs);
    projectile_tree->SetBranchAddress("s0s",&projectile_data.s0s);
    projectile_tree->SetBranchAddress("s0xs",&projectile_data.s0xs);
    projectile_tree->SetBranchAddress("generation",&projectile_data.generation);
    projectile_tree->SetBranchAddress("height",&projectile_data.height);
    projectile_tree->SetBranchAddress("weight",&projectile_data.weight);
    projectile_tree->SetBranchAddress("slantToImpact",&projectile_data.slantToImpact);
    projectile_tree->SetBranchAddress("slantTraversed",&projectile_data.slantTraversed);
    projectile_tree->SetBranchAddress("xShower",&projectile_data.xShower);
    projectile_tree->SetBranchAddress("yShower",&projectile_data.yShower);

    // * particle tree
    particle::data_t particle_data;

    particle_tree->SetBranchAddress("Energy",&particle_data.Energy);
    particle_tree->SetBranchAddress("Px",&particle_data.Px);
    particle_tree->SetBranchAddress("Py",&particle_data.Py);
    particle_tree->SetBranchAddress("Pz",&particle_data.Pz);
    particle_tree->SetBranchAddress("mass",&particle_data.mass);
    particle_tree->SetBranchAddress("x",&particle_data.x);
    particle_tree->SetBranchAddress("y",&particle_data.y);
    particle_tree->SetBranchAddress("time",&particle_data.time);
    particle_tree->SetBranchAddress("id",&particle_data.id);
    particle_tree->SetBranchAddress("interactionCounter",&particle_data.interactionCounter);
    particle_tree->SetBranchAddress("z",&particle_data.z);
    particle_tree->SetBranchAddress("t_formation",&particle_data.t_formation);
    particle_tree->SetBranchAddress("t_destruction",&particle_data.t_destruction);
    particle_tree->SetBranchAddress("id_origin",&particle_data.id_origin);
    particle_tree->SetBranchAddress("id_father",&particle_data.id_father);
    particle_tree->SetBranchAddress("id_mother",&particle_data.id_mother);
    particle_tree->SetBranchAddress("status",&particle_data.status);

    // * get energy of the first interaction and determine threshold energy
    interaction_tree->GetEntry(0);
    const units::energy_t eprim = interaction_data.eProj * 1_GeV;
    const units::energy_t ethreshold = threshold*eprim;

    // * put all data of the current event into memory
    long long ipart = 0;

    m_interactions.reserve(interaction_tree->GetEntries());
    
    for (long long i = 0; i < interaction_tree->GetEntries(); ++i) {
      // * interaction
      // read current interaction + seeds + projectile into the interaction_reader
      interaction_tree->GetEntry(i);
      seed_tree->GetEntry(i);

      // check if projectile energy is below defined fraction of the primary energy
      if (interaction_data.eProj * 1_GeV < ethreshold) {
        // skip secondaries of this interaction
        ipart += interaction_data.mult;
        // skip the current interaction
        continue;
      }

      // create a new interaction from the last read interaction data and add it to
      // the interaction vetor
      interaction_ptr current_interaction = std::make_shared<interaction>(interaction_data);
      m_interactions.emplace_back(current_interaction);

      // * projectile
      // read the current projectile
      projectile_tree->GetEntry(i);

      // force interaction id of the primary particle to be = 0
      if (i == 0) {
        projectile_data.generation = 0;
      }

      // add projectile to the current interaction
      current_interaction->set_projectile(projectile_data);

      // - consistency check: ensure projectile and interaction have the same counter
      if (check && projectile_data.interactionCounter != interaction_data.interactionCounter) {
        throw std::runtime_error("conex extensions projectile/interaction counter mismatch");
      }

      // * secondary particles
      // get iterator for the particles in the current interaction, then loop
      // over all particles of this interaction and add them to the interaction
      // object
      const auto multiplicity = current_interaction->get_multiplicity();
      const long long ipart_end = ipart + multiplicity;

      units::energy_t esum(0);
      util::vector_t<units::momentum_t> psum(current_interaction->get_frame());

      while (ipart < ipart_end) {
        particle_tree->GetEntry(ipart++);

        // - consistency check: ensure partcile and interaction have the same counter
        if (check && particle_data.interactionCounter != interaction_data.interactionCounter) {
          throw std::runtime_error("conex extensions particle/interaction counter mismatch");
        }

        particle_ptr current_particle = current_interaction->add_particle(particle_data);

        // total energy and momentum for consistency check
        if (check) {
          esum += current_particle->get_energy();
          psum += current_particle->get_momentum();
        }
      }

      // - consistency check: sum of secondary energies must match initial lab energy (eProd)
      if (check && std::fabs(esum/current_interaction->get_lab_energy() - 1) > 1e-7 /* fp_precision */) {
        throw std::runtime_error("conex extensions energy mismatch");
      }

      // - consistency check: sum o secondary momenta must match projectile momentum
      if (check) {
        const util::vector_t<units::momentum_t> pproj = current_interaction->get_projectile()->get_momentum();
        const double pdev = (psum-pproj).norm() / pproj.norm();
        if (pdev > 1e-2 /* 1% */) {
          std::stringstream str;
          str << "conex extensions momentum mismatch" << std::endl;
          str << "psum: " << psum.on_frame(util::frame::conex_observer) << std::endl;
          str << "pproj: " << pproj.on_frame(util::frame::conex_observer) << std::endl;
          str << "dev: " << pdev << std::endl;
          throw std::runtime_error(str.str());
        }
      }

    }

    // delete the trees
    delete interaction_tree;
    delete seed_tree;
    delete projectile_tree;
    delete particle_tree;
  }

} // namespace conex::extensions