#include "cherenkov_photon.h"

namespace corsika {

  cherenkov_photon::cherenkov_photon(
    const float* data, const bool has_thinning
  )
  : particle(data, has_thinning),
    n(m_data[0]),
    x(m_data[1]),
    y(m_data[2]),
    u(m_data[3]),
    v(m_data[4]),
    t(m_data[5]),
    h(m_data[5]),
    w(m_data[6])
  {

  }

  cherenkov_photon::cherenkov_photon(
    const particle& part
  )
  : particle(part),
    n(m_data[0]),
    x(m_data[1]),
    y(m_data[2]),
    u(m_data[3]),
    v(m_data[4]),
    t(m_data[5]),
    h(m_data[5]),
    w(m_data[6])
  {
    
  }

  cherenkov_photon::cherenkov_photon(
    particle&& part
  )
  : particle(part),
    n(m_data[0]),
    x(m_data[1]),
    y(m_data[2]),
    u(m_data[3]),
    v(m_data[4]),
    t(m_data[5]),
    h(m_data[5]),
    w(m_data[6])
  {

  }

} // namespace corsika