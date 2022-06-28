#ifndef _conex_extensions_file_h
#define _conex_extensions_file_h

#include <string>
#include <cstddef>

#include <TFile.h>
#include <TTree.h>

#include <conex/extensions/event.h>

namespace conex::extensions {

  class file : private TFile {
  private:
    size_t m_event_count;

  public:
    file(const std::string& fname);

    event get_event(const size_t pos, const double threshold = 0.01, const bool check = false);

    bool is_open() const
    {return IsOpen() && !IsZombie() && !TestBit(kRecovered);}

    size_t get_n_events() const
    {return m_event_count;}
  };

} // namespace conex::extensions

#endif // _conex_extensions_file_h