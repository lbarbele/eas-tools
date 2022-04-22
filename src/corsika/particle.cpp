#include <algorithm>

#include <corsika/particle.h>

namespace corsika {

  particle::particle()
  {
    std::fill(m_data.begin(), m_data.end(), 0.0);
  }

  particle::particle(
    const float* data,
    const bool has_thinning
  )
  {
    if (has_thinning) {
      std::copy(data, data+8, m_data.data());
    } else {
      std::copy(data, data+7, m_data.data());
      m_data.back() = 1;
    }
  }

} // namespace corsika