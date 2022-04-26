#ifndef _conex_file_h
#define _conex_file_h

#include <string>
#include <memory>
#include <cstddef>

#include <TFile.h>
#include <TTree.h>

#include <conex/header.h>
#include <conex/shower.h>

namespace conex {

  class file : private TFile {
  public:
    class iterator;
  // TODO: add support to the leading interactions tree
  // TODO: add a shower iterator, so that range-based for loops can be used
  private:

    std::unique_ptr<TTree> m_header_tree;
    std::unique_ptr<TTree> m_shower_tree;

    header m_header;
    shower m_shower;

  public:
    file(const std::string& fname);

    bool is_open() const;
    const header& get_header() const;
    long long get_n_showers() const;
    const shower& get_shower(size_t pos);
    const shower& get_shower() const;

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