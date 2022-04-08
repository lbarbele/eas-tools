#ifndef _corsika_shower_h
#define _corsika_shower_h

#include "subblock.h"
#include "fstream-iterator.h"

namespace corsika {

  class shower {
  private:
    std::shared_ptr<fstream> m_stream;
    fstream::iterator m_header;
    fstream::iterator m_trailer;
  public:
    shower(const std::shared_ptr<fstream>& stream, fstream::iterator& it);

    const subblock& get_header() const
    {return *m_header;}

    const subblock& get_trailer() const
    {return *m_trailer;}
  };

} // namespace corsika

#endif