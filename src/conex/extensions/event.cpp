#include "event.h"

namespace conex::extensions {

  event::event(
    TTree* particle_tree,
    TTree* projectile_tree,
    TTree* interaction_tree,
    TTree* seed_tree
  ) :
    m_particle_tree(particle_tree),
    m_projectile_tree(projectile_tree),
    m_interaction_tree(interaction_tree),
    m_seed_tree(seed_tree)
  {
    // particle tree
    m_particle_tree->SetBranchAddress("Energy",&m_particle.Energy);
    m_particle_tree->SetBranchAddress("Px",&m_particle.Px);
    m_particle_tree->SetBranchAddress("Py",&m_particle.Py);
    m_particle_tree->SetBranchAddress("Pz",&m_particle.Pz);
    m_particle_tree->SetBranchAddress("mass",&m_particle.mass);
    m_particle_tree->SetBranchAddress("x",&m_particle.x);
    m_particle_tree->SetBranchAddress("y",&m_particle.y);
    m_particle_tree->SetBranchAddress("time",&m_particle.time);
    m_particle_tree->SetBranchAddress("id",&m_particle.id);
    m_particle_tree->SetBranchAddress("interactionCounter",&m_particle.interactionCounter);
    m_particle_tree->SetBranchAddress("z",&m_particle.z);
    m_particle_tree->SetBranchAddress("t_formation",&m_particle.t_formation);
    m_particle_tree->SetBranchAddress("t_destruction",&m_particle.t_destruction);
    m_particle_tree->SetBranchAddress("id_origin",&m_particle.id_origin);
    m_particle_tree->SetBranchAddress("id_father",&m_particle.id_father);
    m_particle_tree->SetBranchAddress("id_mother",&m_particle.id_mother);
    m_particle_tree->SetBranchAddress("status",&m_particle.status);

    // projectile tree
    m_projectile_tree->SetBranchAddress("Energy",&m_interaction.m_projectile.Energy);
    m_projectile_tree->SetBranchAddress("Px",&m_interaction.m_projectile.Px);
    m_projectile_tree->SetBranchAddress("Py",&m_interaction.m_projectile.Py);
    m_projectile_tree->SetBranchAddress("Pz",&m_interaction.m_projectile.Pz);
    m_projectile_tree->SetBranchAddress("mass",&m_interaction.m_projectile.mass);
    m_projectile_tree->SetBranchAddress("x",&m_interaction.m_projectile.x);
    m_projectile_tree->SetBranchAddress("y",&m_interaction.m_projectile.y);
    m_projectile_tree->SetBranchAddress("time",&m_interaction.m_projectile.time);
    m_projectile_tree->SetBranchAddress("id",&m_interaction.m_projectile.id);
    m_projectile_tree->SetBranchAddress("interactionCounter",&m_interaction.m_projectile.interactionCounter);
    m_projectile_tree->SetBranchAddress("c0s",&m_interaction.m_projectile.c0s);
    m_projectile_tree->SetBranchAddress("c0xs",&m_interaction.m_projectile.c0xs);
    m_projectile_tree->SetBranchAddress("s0s",&m_interaction.m_projectile.s0s);
    m_projectile_tree->SetBranchAddress("s0xs",&m_interaction.m_projectile.s0xs);
    m_projectile_tree->SetBranchAddress("generation",&m_interaction.m_projectile.generation);
    m_projectile_tree->SetBranchAddress("height",&m_interaction.m_projectile.height);
    m_projectile_tree->SetBranchAddress("weight",&m_interaction.m_projectile.weight);
    m_projectile_tree->SetBranchAddress("slantToImpact",&m_interaction.m_projectile.slantToImpact);
    m_projectile_tree->SetBranchAddress("slantTraversed",&m_interaction.m_projectile.slantTraversed);
    m_projectile_tree->SetBranchAddress("xShower",&m_interaction.m_projectile.xShower);
    m_projectile_tree->SetBranchAddress("yShower",&m_interaction.m_projectile.yShower);

    // interaction tree
    m_interaction_tree->SetBranchAddress("idProj",&m_interaction.idProj);
    m_interaction_tree->SetBranchAddress("idTarg",&m_interaction.idTarg);
    m_interaction_tree->SetBranchAddress("mult",&m_interaction.mult);
    m_interaction_tree->SetBranchAddress("eProj",&m_interaction.eProj);
    m_interaction_tree->SetBranchAddress("eCMS",&m_interaction.eCMS);
    m_interaction_tree->SetBranchAddress("eProd",&m_interaction.eProd);
    m_interaction_tree->SetBranchAddress("interactionCounter",&m_interaction.interactionCounter);

    // seed tree
    m_seed_tree->SetBranchAddress("seed1a",&m_interaction.seed1);
    m_seed_tree->SetBranchAddress("seed2a",&m_interaction.seed2);
    m_seed_tree->SetBranchAddress("seed3a",&m_interaction.seed3);

    // get indices of the particle tree
    m_particle_tree_indices.push_back(0);
    for (int i = 0; i < get_n_interactions(); ++i) {
      m_interaction_tree->GetEntry(i);
      m_particle_tree_indices.push_back(
        m_particle_tree_indices.back() + m_interaction.mult
      );
    }

  }

  long long
  event::get_n_interactions()
  const
  {
    return m_interaction_tree->GetEntries();
  }

  interaction
  event::get_interaction(
    size_t pos
  )
  {
    m_projectile_tree->GetEntry(pos);
    m_interaction_tree->GetEntry(pos);
    m_seed_tree->GetEntry(pos);

    interaction ret = m_interaction;

    auto index = m_particle_tree_indices[pos];
    auto past_the_end = m_particle_tree_indices[pos+1];

    while (index < past_the_end) {
      m_particle_tree->GetEntry(index++);
      ret.m_secondaries.push_back(m_particle);
    }

    return ret;
  }

} // namespace conex::extensions