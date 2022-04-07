
#include "shower.h"
#include "fstream-iterator.h"

namespace corsika {

  shower::shower(
    const fstream::iterator& it
  ) :
    m_begin(it)
  {
  }

} // namespace corsika