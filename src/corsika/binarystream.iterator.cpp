#include <corsika/binarystream.h>
#include <corsika/subblock.h>

#include <iostream>

namespace corsika {

  binarystream::iterator::iterator()
  : m_stream(nullptr),
    m_pos(-1)
  {}

  binarystream::iterator::iterator(
    binarystream& stream
  ) :
    m_stream(&stream)
  {
    stream.seekg(0);

    switch (stream.get<int32_t>())
    {
    case 22932:
      m_subblock_size = 273;
      m_block_size = 22932/sizeof(binarystream::char_type);
      m_has_thinning = false;
      break;
    case 26208:
      m_subblock_size = 312;
      m_block_size = 26208/sizeof(binarystream::char_type) + 2;
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
  binarystream::iterator::get_pos()
  const
  {
    return m_pos;
  }

  bool
  binarystream::iterator::has_thinning()
  const
  {
    return m_has_thinning;
  }

  const subblock& 
  binarystream::iterator::operator*()
  const
  {
    return m_data;
  }

  const subblock*
  binarystream::iterator::operator->()
  const
  {
    return &m_data;
  }

  binarystream::iterator&
  binarystream::iterator::operator++()
  {
    ++m_pos;

    // position the stream should be
    const long pos_it = (m_pos/21)*m_block_size + 1 + (m_pos%21)*m_subblock_size;
    // position the stream actually is
    const long pos_str = m_stream->tellg()/sizeof(binarystream::char_type);

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

  binarystream::iterator
  binarystream::iterator::operator++(int)
  {
    binarystream::iterator other(*this);
    ++(*this);
    return other;
  }

  binarystream::iterator
  binarystream::iterator::operator+(
    long offset
  ) const
  {
    if (offset < 0) {
      std::cerr << "WARNING: decrement is not supported by corsika::fstream::iterator!" << std::endl;
      throw;
    }

    auto end = binarystream::iterator();
    auto ret = *this;

    while (offset > 0 && ret != end) {
      ++ret;
      --offset;
    }

    return ret;
  }

  bool 
  binarystream::iterator::operator==(
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
  binarystream::iterator::operator!=(
    const iterator& other
  ) const
  {
    return !this->operator==(other);
  }

} // namespace corsika