#ifndef _corsika_shower_h
#define _corsika_shower_h

#include <corsika/subblock.h>
#include <corsika/binarystream.h>
#include <corsika/particle_iterator.h>

namespace corsika {

  class shower {
  private:
    std::shared_ptr<binarystream> m_stream;
    binarystream::iterator m_header;
    binarystream::iterator m_trailer;
  public:
    shower(const std::shared_ptr<binarystream>& stream, binarystream::iterator& it);

    const subblock& get_header() const
    {return *m_header;}

    const subblock& get_trailer() const
    {return *m_trailer;}

    particle_iterator begin() const
    {return particle_iterator(m_header+1);}

    particle_iterator end() const
    {return particle_iterator();}

  };

} // namespace corsika

#endif