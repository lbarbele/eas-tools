#ifndef _corsika_binaryfile_h
#define _corsika_binaryfile_h

#include <corsika/shower.h>
#include <corsika/binarystream.h>
#include <corsika/subblock.h>

#include <vector>
#include <memory>
#include <string>

namespace corsika {

  class binaryfile {
  private:
    std::shared_ptr<binarystream> m_stream;
    subblock m_header;
    subblock m_trailer;
    std::vector<shower> m_showers;

  public:
    binaryfile(const std::string& fname);

    const subblock& get_header() const
    {return m_header;}

    const subblock& get_trailer() const
    {return m_trailer;}

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