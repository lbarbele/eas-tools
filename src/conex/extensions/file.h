#ifndef _conex_extensions_particlefile_h
#define _conex_extensions_particlefile_h

#include <string>
#include <cstddef>

#include <TFile.h>
#include <TTree.h>

#include "event.h"
#include "particle.h"
#include "projectile.h"

namespace conex::extensions {

  class file : private TFile {
  private:
    size_t m_event_count;

  public:
    file(const std::string& fname);

    bool is_open() const;
    event get_event(size_t pos);
    size_t get_n_events() const;
  };

} // namespace conex::extensions

#endif // _conex_particle_file_h