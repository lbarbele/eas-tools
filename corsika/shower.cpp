
#include <iostream>

#include "shower.h"
#include "fstream-iterator.h"

namespace corsika {

  shower::shower(
    const std::shared_ptr<fstream>& stream,
    fstream::iterator& it
  ) :
    m_stream(stream),
    m_begin(it)
  {
    // check if iterator is actually poiting to the event header block
    if (it->get_type() != subblock::type::event_header) {
      std::cerr << "corsika::shower::shower(): iterator is not an event header!" << std::endl;
      throw;
    }

    // consume the stream until the event end block is found or the stream becomes unreadable
    // TODO: add support to LONG blocks
    while(it->get_type() != subblock::type::event_end && it != stream->end()) {
      ++it;
    }

    // ensure we got an event end block
    if (it->get_type() == subblock::type::event_end) {
      m_end = it;
    } else {
      std::cerr << "corsika::shower::shower(): unable to find event trailer block!" << std::endl;
      throw;
    }
  }

} // namespace corsika