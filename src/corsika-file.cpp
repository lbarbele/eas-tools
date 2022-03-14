#include <string>
#include <memory>
#include <iostream>

#include "corsika-file.h"
#include "corsika-fstream.h"
#include "corsika-fstream-iterator.h"

namespace corsika {

  file::file(
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
        std::cout << "found event header at " << it.get_pos() << std::endl;
        break;
      case subblock::type::event_end:
        std::cout << "found event end at " << it.get_pos() << std::endl;
        break;
      case subblock::type::longitudinal:
       std::cout << "found longitudinal block at " << it.get_pos() << std::endl;
        break;
      case subblock::type::data:
        break;
      default:
        throw;
      }
    }
  }

} // namespace corsika