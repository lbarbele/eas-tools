#include <corsika/profile.h>

namespace corsika {

  profile::profile()
  : m_is_slant(false),
    m_size(0),
    m_step_size(-1)
  {
  }

  void
  profile::clear()
  {
    for (auto& v : m_data) {
      v.resize(0);
    }

    m_size = 0;
  }

  void
  profile::resize(
    const unsigned int size
  )
  {
    for (auto& v : m_data) {
      v.resize(size, 0);
    }

    m_size = size;
  }

} // namespace corsika