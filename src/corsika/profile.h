#ifndef _corsika_profile_h
#define _corsika_profile_h

#include <util/gaisser_hillas_fit.h>

#include <ostream>
#include <iomanip>
#include <array>
#include <valarray>

namespace corsika {

  class profile {
  public:
    enum class type : unsigned int {
      depth = 0, gammas, positrons, electrons, mu_plus, mu_minus, hadrons, charged, nuclei,
      depth_dep, gamma_dep, em_ioniz, em_cut, mu_ioniz, mu_cut, hadr_ioniz, hadr_cut, netrino, dedx_sum
    };

  private:
    bool m_is_slant;
    unsigned int m_size;
    double m_step_size;
    util::gaisser_hillas_fit m_fit;
    std::array<std::valarray<double>, 20> m_data;

  public:
    profile();

    void clear();
    void resize(const unsigned int size);

    std::valarray<double>& get(const unsigned int iprof)
    {return m_data.at(iprof);}
    const std::valarray<double>& get(const unsigned int iprof) const
    {return m_data.at(iprof);}

    std::valarray<double>& get(const type tp)
    {return get(static_cast<unsigned int>(tp));}
    const std::valarray<double>& get(const type tp) const
    {return get(static_cast<unsigned int>(tp));}

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
          stream << std::setw(13) << prof.get(icol)[istep];
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