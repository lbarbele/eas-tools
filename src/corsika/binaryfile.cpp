#include <string>
#include <memory>
#include <iostream>

#include "binaryfile.h"
#include "fstream.h"
#include "fstream-iterator.h"

namespace corsika {

  binaryfile::binaryfile(
    const std::string& fname
  ) :
    m_stream(std::make_unique<fstream>(fname, std::ios::in))
  {
    for (auto it = m_stream->begin(); it != m_stream->end(); ++it) {
      switch (it->get_type())
      {
      case subblock::type::run_header:
        m_header = *it;
        break;
      case subblock::type::run_end:
        m_end = *it;
        break;
      case subblock::type::event_header:
        m_showers.emplace_back(it);
        return;
      case subblock::type::event_end:
        break;
      case subblock::type::longitudinal:
        break;
      case subblock::type::data:
        break;
      default:
        throw;
      }
    }
  }

} // namespace corsika