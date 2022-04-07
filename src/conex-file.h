#ifndef _conex_file_h
#define _conex_file_h

#include <string>
#include <memory>
#include <cstddef>

#include <TFile.h>
#include <TTree.h>

#include "conex-header.h"
#include "conex-shower.h"

namespace conex {

  class file : private TFile {
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

    const shower& operator[](size_t pos);

  };

} // namespace conex

#endif // _conex_file_h