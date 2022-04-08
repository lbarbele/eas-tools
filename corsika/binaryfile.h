#ifndef _corsika_file_h
#define _corsika_file_h

#include <vector>
#include <memory>
#include <string>

#include "shower.h"
#include "fstream.h"
#include "subblock.h"

namespace corsika {

  class binaryfile {
  private:
    std::shared_ptr<fstream> m_stream;
    subblock m_header;
    subblock m_trailer;
    std::vector<shower> m_showers;

  public:
    binaryfile(const std::string& fname);

    const subblock& get_header() const
    {return m_header;}

    const subblock& get_trailer() const
    {return m_trailer;}

    // TODO: add support to reverse iterators
    auto begin(){return m_showers.begin();}
    auto cbegin() const {return m_showers.cbegin();}
    auto rbegin() {return m_showers.rbegin();}
    auto crbegin() const {return m_showers.crbegin();}
    auto end() {return m_showers.end();}
    auto cend() const {return m_showers.cend();}
    auto rend() {return m_showers.rend();}
    auto crend() const {return m_showers.crend();}
  };

} // namespace corsika 

#endif // _corsika_file_h