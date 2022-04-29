#ifndef _corsika_particle_iterator_h
#define _corsika_particle_iterator_h

#include <corsika/particle.h>
#include <corsika/subblock.h>
#include <corsika/binarystream.h>

namespace corsika {

  class particle_iterator {
  public:
    binarystream::iterator m_stream_it;
    subblock::const_iterator m_subblock_it;
    particle m_particle;

  public:
    particle_iterator();
    particle_iterator(const binarystream::iterator& it);

    const particle& operator*() const;
    const particle* operator->() const;
    particle_iterator& operator++();
    particle_iterator operator++(int);

    bool operator==(const particle_iterator& other) const;
    bool operator!=(const particle_iterator& other) const;
  };

} // namespace corsika


#endif // _corsika_particle_iterator_h