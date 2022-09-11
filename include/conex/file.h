#ifndef _conex_file_h
#define _conex_file_h

#include <string>
#include <vector>
#include <memory>
#include <cstddef>
#include <istream>
#include <iterator>
#include <concepts>

#include <TChain.h>

#include <conex/header.h>
#include <conex/shower.h>

namespace conex {

  class file {
  public:
    class iterator;
  // TODO: add support to the leading interactions tree
  private:
    std::unique_ptr<TChain> m_shower_tree;
    shower m_shower;
    std::vector<std::string> m_file_names;

  public:
    // main constructor: receives a vector of file names
    file(const std::vector<std::string>& fnames);

    // constructor receiving a list of (convertible to) string arguments, each one pointing to a
    // single file. delegates construction to the main constructor
    template <std::convertible_to<std::string>... Args>
    file(Args... args)
    : file(std::vector<std::string>{std::string(args)...})
    {}

    // constructor receiving an input stream. reads the stream contents as a sequence of strings
    // into a vector and delegates to the main constructor
    template <class CharT, class Traits>
    file(std::basic_istream<CharT, Traits>& stream)
    : file(std::vector<std::string>{std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>()})
    {}

    bool is_open() const;
    header get_header() const;

    unsigned int get_n_files() const;
    int get_current_file_number() const;
    std::string get_current_file_name() const;

    long long get_n_showers() const;
    const shower& get_shower(size_t pos);
    const shower& operator[](size_t pos);

    iterator begin();
    iterator end();

  };

  class file::iterator {
  private:
    int m_pos;
    file* m_file;

  public:
    iterator();
    iterator(file& f);

    int get_pos() const
    {return m_pos;}

    const shower& operator*();
    const shower* operator->();

    iterator& operator++();
    iterator operator++(int);

    bool operator==(const iterator& other) const;
    bool operator!=(const iterator& other) const;
  
  };

} // namespace conex

#endif // _conex_file_h