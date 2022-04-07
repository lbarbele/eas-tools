#ifndef _corsika_shower_h
#define _corsika_shower_h

#include "corsika-subblock.h"
#include "corsika-fstream-iterator.h"

namespace corsika {

  class shower {
  private:
    fstream::iterator m_begin;
  public:
    shower(const fstream::iterator& it);

    const subblock& get_header()
    {return *m_begin;}
  };

} // namespace corsika

#endif