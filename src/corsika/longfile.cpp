#include <corsika/longfile.h>

#include <iostream>
#include <iomanip>
#include <string>

namespace corsika {

  longfile::longfile(
    const std::string& fname
  ) :
    m_stream(std::make_shared<std::ifstream>(fname))
  {
    // check if the underlying stream is open
    if (!m_stream->is_open()) {
      m_stream = nullptr;
      return;
    }
  }

  bool
  longfile::is_open()
  const
  {
    return m_stream && m_stream->is_open();
  }

  longfile::iterator
  longfile::begin()
  {
    return longfile::iterator(m_stream);
  }

  longfile::iterator
  longfile::end()
  {
    return longfile::iterator();
  }

} // namespace corsika