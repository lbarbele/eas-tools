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
      v.clear();
    }

    m_size = 0;
  }

  void
  profile::resize(
    const unsigned int size
  )
  {
    for (auto& v : m_data) {
      v.resize(size);
    }

    m_size = size;
  }

  std::vector<double>&
  profile::get(
    const unsigned int iprof
  )
  {
    return m_data.at(iprof);
  }

  const std::vector<double>&
  profile::get(
    const unsigned int iprof
  ) const
  {
    return m_data.at(iprof);
  }

  // void
  // profile::set_fit(
  //   util::gaisser_hillas_fit&& fit
  // )
  // {
  //   m_fit = fit;
  // }

  // const util::gaisser_hillas_fit&
  // profile::get_fit()
  // const
  // {
  //   return m_fit;
  // }

} // namespace corsika