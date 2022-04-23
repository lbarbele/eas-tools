#ifndef _corsika_profile_h
#define _corsika_profile_h

#include <util/gaisser_hillas_fit.h>

#include <ostream>
#include <iomanip>
#include <array>
#include <vector>

namespace corsika {

  class profile {
  private:
    bool m_is_slant;
    unsigned int m_size;
    double m_step_size;
    util::gaisser_hillas_fit m_fit;
    std::array<std::vector<double>, 20> m_data;

  public:
    profile();

    void clear();
    void resize(const unsigned int size);

    std::vector<double>& get(const unsigned int iprof);
    const std::vector<double>& get(const unsigned int iprof) const;

    bool empty() const
    {return m_size == 0;}

    unsigned int size() const
    {return m_size;}

    void set_slant(const bool is_slant)
    {m_is_slant = is_slant;}
    bool is_slant() const
    {return m_is_slant;}

    void set_step_size(const double step)
    {m_step_size = step;}
    double get_step_size() const
    {return m_step_size;}

    void set_fit(const util::gaisser_hillas_fit& fit)
    {m_fit = fit;}
    const util::gaisser_hillas_fit& get_fit() const
    {return m_fit;}
  };

  // function to print a profile
  template<class CharT, class Traits>
  std::basic_ostream<CharT, Traits>&
  operator<<(
    std::basic_ostream<CharT, Traits>& stream,
    const profile& prof
  )
  {
    for (unsigned int istart : {0, 10}) {
      for (unsigned int istep = 0; istep < prof.size(); ++istep) {
        for (unsigned int icol = istart; icol < istart+10; ++icol) {
          stream << std::setw(13) << prof.get(icol).at(istep);
        }
        stream << std::endl;
      }

      if (istart == 0) {
        stream << std::endl;
      }
    }
    return stream;
  }

};

#endif