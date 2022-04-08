#include <iostream>
#include <iomanip>
#include <list>
#include <string>

// #include "conex-file.h"
// #include "conex-extensions-particlefile.h"
// #include "conex-extensions-interaction.h"
// #include "utl-vector.h"

#include "corsika/binaryfile.h"
#include "corsika/fstream.h"

int main(void)
{
  const std::string baseFileName = "/home/luan/Projetos/offline-tests/reconstructionTests/input/heatSimFiles/DAT013733";
  const std::string cherFileName = baseFileName + ".cher-tel005";
  const std::string longFileName = baseFileName + ".long";
  const std::string lstFileName = baseFileName + ".lst";

  corsika::binaryfile cherFile(cherFileName);

  auto shower = *cherFile.begin();

  for (auto part : shower) {
    std::cout << part[0] << std::endl;
  }

  // auto first = shower.begin().m_stream_it;
  // int n = 0;

  // while(first->get_type() == corsika::subblock::type::data) {
  //   auto part_it = first->begin();
  //   while(part_it != first->end()) {
  //     std::cout << part_it->as<float>() << std::endl;
  //     ++n;
  //     part_it += 8;
  //   }
  //   ++first;
  // }

  // std::cout << n << std::endl;
  // std::cout << first.get_pos() - 1 << std::endl;

  // auto beg = shower.begin();
  // auto beg_another = beg;
  // auto end = shower.end();
  // auto end_another = end;

  // if (beg == beg_another) {
  //   std::cout << "is equal" << std::endl;
  // } else {
  //   std::cout << "is different" << std::endl;
  // }

  return 0;

  // for (; beg != end; ++beg) {

  // }

  // auto stream_it = beg.m_stream_it;
  
  // while(stream_it->get_type() == corsika::subblock::type::data) {
  //   for (auto sb_it = stream_it->begin(); sb_it != stream_it->end(); sb_it+=8) {
  //     corsika::particle part(&sb_it->as<float>());
  //     std::cout << part[0] << std::endl;
  //   }
  //   ++stream_it;
  // }

  // for (auto shower : cherFile) {
  //   auto beg = shower.begin();
  //   auto end = shower.end();
  //   auto it = beg;
  //   for (; it != end; ++it) { auto particle = *it;
  //   // for (auto particle : shower) {
  //   // for (; beg != end; ++beg) { auto particle = *beg;
  //     for (int i = 0; i < 8; ++i) {
  //       std::cout << std::setw(13) << particle[i];
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  return 0;
}

// using utl::vector;

// bool
// is_same(
//   const conex::extensions::particle& part,
//   const conex::extensions::projectile& proj,
//   const double phi,
//   const double theta
// )
// {
//   // check for same ID
//   if (part.get_id() != proj.get_id()) {
//     return false;
//   }

//   // check for energy
//   if (std::fabs(part.get_energy()/proj.get_energy() - 1) > 0.001) {
//     return false;
//   }

//   // check momentum direction
//   auto p_part = part.get_momentum().to_obs(phi, theta);
//   auto p_proj = proj.get_momentum();
//   if ((p_part*p_proj)/(p_part.norm()*p_proj.norm()) < 0.999) {
//     return false;
//   }

//   return true;
// }

// struct stackin_particle {

//   int id;
//   double energy;
//   double px;
//   double py;
//   double pz;
//   double x;
//   double y;
//   double h;
//   double t;
//   int gen;
//   double h_previous;

// };

// int
// main(
//   int argc,
//   char const** argv
// )
// {
//   const double pi = std::acos(-1);
//   const double c = 299792458.;

//   conex::file cxFile("../../double-bump/data/anom_sibyll23d_997649511_100.root");
//   auto shower = cxFile.get_shower(0);
//   double theta = shower.zenith*pi/180.;
//   double phi = shower.azimuth*pi/180.;

//   conex::extensions::particlefile partFile("../../double-bump/data/anom_sibyll23d_997649511_100_part.root");

//   auto event = partFile.get_event(0);

//   // get the first interaction and the primary projectile
//   auto firstInteraction = event.get_interaction(0);
//   conex::extensions::projectile primaryParticle = firstInteraction.get_projectile();

//   // set an energy cutoff above which we will store particles
//   const double energyCut = 0.005 * primaryParticle.Energy;

//   std::list<conex::extensions::particle> particleStack;
//   std::list<conex::extensions::particle> leadingStack;

//   // get the particles from the first interaction
//   for (auto& particle : firstInteraction) {
//     if (particle.Energy >= energyCut) {
//       leadingStack.push_back(particle);
//     } else {
//       particleStack.push_back(particle);
//     }
//   }

//   std::cout << "looking for the interactions of " << leadingStack.size() << " leading particles" << std::endl;

//   // now search for the interactions of these leading particles
//   for (auto leading : leadingStack) {

//     std::cout
//       << "looking for id " << leading.get_id()
//       << " with energy " << leading.get_energy()
//       << " GeV ... ";

//     bool found = false;

//     for (int i = 1; i < event.get_n_interactions(); ++i) {
//       auto interaction = event.get_interaction(i);
//       auto projectile = interaction.get_projectile();

//       if (is_same(leading, projectile, phi, theta)) {
//         std::cout << "found at interaction " << interaction.get_interaction_counter() << std::endl;
//         found = true;
//         break;
//       }
//     }

//     if (!found) {
//       std::cout << "not found!" << std::endl;
//     }
//   }

//   return 0;
// }
