#ifndef _corsika_longfile_h
#define _corsika_longfile_h

#include <corsika/profile.h>

#include <string>
#include <fstream>
#include <memory>
#include <list>

namespace corsika {

  class longfile {
  public:
    class iterator;

  private:
    std::shared_ptr<std::ifstream> m_stream;

  public:
    longfile(const std::string& fname);
    bool is_open() const;

    iterator begin();
    iterator end();
  };

  class longfile::iterator {
  private:
    int m_ishower;
    std::shared_ptr<std::ifstream> m_stream;
    profile m_profile;

    bool read_block();

  public:
    iterator();
    iterator(const std::shared_ptr<std::ifstream>& stream);

    int get_ishower() const
    {return m_ishower;}

    const profile& operator*() const;
    const profile* operator->() const;
    iterator& operator++();
    iterator operator++(int);
    
    bool operator==(const iterator& other) const;
    bool operator!=(const iterator& other) const;
  };

} // namespace corsika

#endif