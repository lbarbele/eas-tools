#include <iostream>

#include "corsika-fstream.h"
#include "corsika-fstream-iterator.h"

namespace corsika {

  fstream::iterator::iterator()
  : m_stream(nullptr)
  {}

  fstream::iterator::iterator(
    fstream& stream
  ) :
    m_stream(&stream)
  {
    stream.seekg(0);
    auto block_size = stream.get<int32_t>();
    switch (block_size)
    {
    case 22932:
      m_subblock_size = block_size / (4 * 21);
      break;
    case 26208:
      m_subblock_size = block_size / (4 * 21);
      break;
    default:
      throw;
    }

    m_pos = 0;
    m_data.resize(m_subblock_size);
    stream.read(m_data.data(), m_subblock_size);
  }

  long
  fstream::iterator::get_pos()
  const
  {
    return m_pos;
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
    if ((++m_pos)%21 == 0) {
      m_stream->ignore(2);
    }

    if (m_stream->good()) {
      m_data.resize(m_subblock_size);
      m_stream->read(m_data.data(), m_subblock_size);
    }

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