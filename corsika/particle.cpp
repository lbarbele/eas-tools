#include <algorithm>

#include "particle.h"

namespace corsika {

  particle::particle()
  {
    std::fill(m_data.begin(), m_data.end(), 0.0);
  }

  particle::particle(
    const float* data
  )
  {
    for (int i = 0; i < 8; i++) {
      m_data[i] = data[i];
    }
  }

} // namespace corsika