#include <iostream>

#include "fstream.h"
#include "fstream-iterator.h"
#include "subblock.h"

namespace corsika {

  fstream::iterator::iterator()
  : m_stream(nullptr),
    m_pos(-1)
  {}

  fstream::iterator::iterator(
    fstream& stream
  ) :
    m_stream(&stream)
  {
    stream.seekg(0);

    switch (stream.get<int32_t>())
    {
    case 22932:
      m_subblock_size = 273;
      m_block_size = 22932/sizeof(fstream::char_type);
      m_has_thinning = false;
      break;
    case 26208:
      m_subblock_size = 312;
      m_block_size = 26208/sizeof(fstream::char_type) + 2;
      m_has_thinning = true;
      break;
    default:
      throw;
    }

    m_pos = 0;
    m_data.resize(m_subblock_size);
    stream.read(m_data.data(), m_subblock_size);

    if (!m_stream->good() || m_stream->gcount() != m_subblock_size) {
      m_data.clear();
      m_stream = nullptr;
      m_pos = -1;
    }
  }

  long
  fstream::iterator::get_pos()
  const
  {
    return m_pos;
  }

  bool
  fstream::iterator::has_thinning()
  const
  {
    return m_has_thinning;
  }

  const subblock& 
  fstream::iterator::operator*()
  const
  {
    return m_data;
  }

  const subblock*
  fstream::iterator::operator->()
  const
  {
    return &m_data;
  }

  fstream::iterator&
  fstream::iterator::operator++()
  {
    ++m_pos;

    // position the stream should be
    const long pos_it = (m_pos/21)*m_block_size + 1 + (m_pos%21)*m_subblock_size;
    // position the stream actually is
    const long pos_str = m_stream->tellg()/sizeof(fstream::char_type);

    // reposition the stream if necessary
    if (pos_it != pos_str) {
      m_stream->clear();
      m_stream->seekg(pos_it, std::ios::beg);
    }

    // only perform the actual reading if the stream is still good
    if (m_stream->good()) {
      m_stream->read(m_data.data(), m_subblock_size);
    }

    // if the stream is not good anymore, let this iterator be the end of
    // stream iterator
    if (!m_stream->good() || m_stream->gcount() != m_subblock_size) {
      m_data.clear();
      m_stream = nullptr;
      m_pos = -1;
    }
    
    return *this;
  }

  fstream::iterator
  fstream::iterator::operator++(int)
  {
    fstream::iterator other(*this);
    ++(*this);
    return other;
  }

  fstream::iterator
  fstream::iterator::operator+(
    long offset
  ) const
  {
    if (offset < 0) {
      std::cerr << "WARNING: decrement is not supported by corsika::fstream::iterator!" << std::endl;
      throw;
    }

    auto end = fstream::iterator();
    auto ret = *this;

    while (offset > 0 && ret != end) {
      ++ret;
      --offset;
    }

    return ret;
  }

  bool 
  fstream::iterator::operator==(
    const iterator& other
  ) const
  {
    return !m_stream?
      other.m_stream == nullptr : (
        m_stream == other.m_stream &&
        m_pos == other.m_pos
      );
  }

  bool
  fstream::iterator::operator!=(
    const iterator& other
  ) const
  {
    return !this->operator==(other);
  }

} // namespace corsika