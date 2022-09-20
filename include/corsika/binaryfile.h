#ifndef _corsika_binaryfile_h
#define _corsika_binaryfile_h

#include <corsika/shower.h>
#include <corsika/binarystream.h>
#include <corsika/subblock.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace corsika {

  class binaryfile {
  private:
    std::shared_ptr<binarystream> m_stream;
    subblock m_header;
    subblock m_trailer;
    std::vector<shower> m_showers;

  public:
    binaryfile(const std::string& fname)
    : m_stream(std::make_shared<binarystream>(fname, std::ios::in))
    {
      // consume the stream until the event end is found
      for (auto it = m_stream->begin(); it != m_stream->end(); ++it) {
        switch (it->get_type())
        {
        case subblock::type::run_header:
          m_header = *it;
          continue;
        case subblock::type::run_end:
          m_trailer = *it;
          return;
        case subblock::type::event_header:
          m_showers.emplace_back(m_stream, it);
          continue;
        default:
          std::cerr << "corsika::binaryfile::binaryfile(): unexpected block " << it->get_title() << std::endl;
          throw;
        }
      }

      // if we get here, the event trailer was not found
      std::cerr << "corsika::binaryfile::binaryfile(): could not found the event trailer block!" << std::endl;
      throw;
    }

    const subblock& get_header() const
    {return m_header;}

    const subblock& get_trailer() const
    {return m_trailer;}

    auto get_n_showers() const
    {return m_showers.size();}

    const shower& get_shower(const std::size_t i) const
    {return m_showers.at(i);}

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