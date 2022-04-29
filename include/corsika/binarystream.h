#ifndef _corsika_binarystream_h
#define _corsika_binarystream_h

#include <corsika/word_t.h>
#include <corsika/subblock.h>

#include <fstream>
#include <filesystem>
#include <locale>
#include <string>

namespace corsika {

  class binarystream : public std::basic_fstream<word_t> {
  private:
  public:
    class iterator;

    binarystream();
    binarystream(const char* fname, ios_base::openmode mode = ios_base::in);
    binarystream(const std::string& fname, ios_base::openmode mode = ios_base::in);
    binarystream(const std::filesystem::path& fname, ios_base::openmode mode = ios_base::in);

    void open(const char* fname, ios_base::openmode mode = ios_base::in);
    void open(const std::string& fname, ios_base::openmode mode = ios_base::in);
    void open(const std::filesystem::path& fname, ios_base::openmode mode = ios_base::in);

    iterator begin();
    iterator end();

    template<typename T = word_t>
    T get()
    {
      word_t word;
      std::basic_fstream<word_t>::get(word);
      return word.as<T>();
    }
  };

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

  // word insertion/extraction operators
  std::basic_istream<word_t>&
  operator>>(std::basic_istream<word_t>& stream, word_t& word);

  std::basic_ostream<word_t>&
  operator<<(std::basic_ostream<word_t>& stream, word_t& word);

} // namespace corsika

#endif