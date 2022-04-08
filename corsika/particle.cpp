#include <algorithm>

#include "particle.h"

namespace corsika {

  particle::particle(
    const float* data
  )
  {
    std::copy(m_data.begin(), m_data.end(), data);
  }

} // namespace corsika