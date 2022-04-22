#ifndef _corsika_fstream_iterator_h
#define _corsika_fstream_iterator_h

#include <corsika/binarystream.h>
#include <corsika/subblock.h>

#include <memory>
#include <vector>

namespace corsika {
  
  class binarystream::iterator {
  private:
    binarystream* m_stream;
    subblock m_data;
    long m_block_size;
    long m_subblock_size;
    long m_pos;
    bool m_has_thinning;

  public:
    iterator();
    iterator(binarystream& stream);

    long get_pos() const;
    bool has_thinning() const;

    const subblock& operator*() const;
    const subblock* operator->() const;
    iterator& operator++();
    iterator operator++(int);
    iterator operator+(long offset) const;

    bool operator==(const iterator& other) const;
    bool operator!=(const iterator& other) const;
  };

} // namespace corsika

#endif