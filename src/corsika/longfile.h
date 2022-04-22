#ifndef _corsika_longfile_h
#define _corsika_longfile_h

#include <string>
#include <fstream>
#include <memory>

namespace corsika {

  // enum class profile_type : int {
  //   depth, gamma, positron, electron, muplus, muminus, hadron, charged, nuclei, cherenkov, all
  // };

  class longfile {
  private:
    std::shared_ptr<std::ifstream> m_stream;
    // std::map<profile_type, std::vector<double>> m_particle_profiles;
    // std::map<profile_type, std::vector<double>> m_deposit_profiles;
    // gaisser_hillas_fit m_fit;

    // static std::vector<profile_type> get_ordered_types()
    // {
    //   // in the same order they appear in CORSIKA!
    //   using t = profile_type;
    //   return {t::depth, t::gamma, t::positron, t::electron, t::muplus, t::muminus, t::hadron, t::charged, t::nuclei, t::cherenkov};
    // }

  public:
    longfile(const std::string& fname);
    bool is_open() const;

    bool read();

    // const std::vector<double>& get_particle_profile(profile_type type) const;
    // const std::vector<double>& get_deposit_profiles(profile_type type) const;

    // const gaisser_hillas_fit& get_fit() const
    // {return m_fit;}

  };

} // namespace corsika

#endif // _corsika_longfile_h