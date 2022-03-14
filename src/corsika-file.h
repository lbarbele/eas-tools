#ifndef _corsika_file_h
#define _corsika_file_h

#include <memory>
#include <string>

#include "corsika-fstream.h"
#include "corsika-subblock.h"

namespace corsika {

  class file {
  private:
    std::shared_ptr<fstream> m_stream;
    subblock m_header;
    subblock m_end;

  public:
    file(const std::string& fname);

    const subblock& get_header() const
    {return m_header;}

    const subblock& get_end() const
    {return m_end;}
  };

} // namespace corsika 

#endif // _corsika_file_h