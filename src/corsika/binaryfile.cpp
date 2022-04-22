#include <string>
#include <memory>
#include <iostream>

#include <corsika/binaryfile.h>
#include <corsika/fstream.h>
#include <corsika/fstream-iterator.h>

namespace corsika {

  binaryfile::binaryfile(
    const std::string& fname
  ) :
    m_stream(std::make_shared<fstream>(fname, std::ios::in))
  {
    // consume the stream until the event end is found
    for (auto it = m_stream->begin(); it != m_stream->end(); ++it) {
      switch (it->get_type())
      {
      case subblock::type::run_header:
        m_header = *it;
        continue;
      case subblock::type::run_end:
        m_trailer = *it;
        return;
      case subblock::type::event_header:
        m_showers.emplace_back(m_stream, it);
        continue;
      default:
        std::cerr << "corsika::binaryfile::binaryfile(): unexpected block " << it->get_title() << std::endl;
        throw;
      }
    }

    // if we get here, the event trailer was not found
    std::cerr << "corsika::binaryfile::binaryfile(): could not found the event trailer block!" << std::endl;
    throw;
  }

} // namespace corsika