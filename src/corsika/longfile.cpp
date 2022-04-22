#include <iostream>
#include <limits>
#include <iomanip>

#include <corsika/longfile.h>

namespace corsika {

  longfile::longfile(
    const std::string& fname
  ) :
    m_stream(std::make_shared<std::ifstream>(fname))
  {
    if (!m_stream->is_open()) {
      m_stream = nullptr;
      return;
    }
  }

  bool
  longfile::is_open()
  const
  {
    return m_stream && m_stream->is_open();
  }

  // bool
  // longfile::read()
  // {
  //   // get a reference to the stream (easier to read)
  //   auto& stream = *m_stream;

  //   //
  //   // particle profiles
  //   // 

  //   // get number of steps in the profiles
  //   unsigned int steps;
  //   stream.ignore(29);
  //   stream >> steps;

  //   // (re)allocate the containers
  //   for (auto type : get_ordered_types()) {
  //     m_particle_profiles[type] = std::vector<double>(steps, 0);
  //   }

  //   // ignore rest of the first line
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //   // ignore line with headers
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // parse the longitudinal particle profiles
  //   for (unsigned int i = 0; i < steps; ++i) {
  //     for (auto type : get_ordered_types()) {
  //       stream >> m_particle_profiles[type][i];
  //     }
  //   }

  //   //
  //   // energy deposit profiles
  //   //

  //   // consume up to the EOL char
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // read number of steps again
  //   stream.ignore(31);
  //   stream >> steps;

  //   // (re)allocate the containers
  //   for (auto type : get_ordered_types()) {
  //     m_deposit_profiles[type] = std::vector<double>(steps, 0);
  //   }

  //   // ignore rest of the first line
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //   // ignore line with headers
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // parse the longitudinal energy deposit profiles
  //   for (unsigned int i = 0; i < steps; ++i) {
  //     for (auto type : get_ordered_types()) {
  //       stream >> m_deposit_profiles[type][i];
  //     }
  //   }

  //   //
  //   // gaisser hillas parameters
  //   //
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  //   stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // stream.ignore(21);
  //   // stream >> m_fit.nmax >> m_fit.x0 >> m_fit.xmax >> m_fit.p1 >> m_fit.p2 >> m_fit.p3;
  //   // stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // stream.ignore(21);
  //   // stream >> m_fit.chi2perdof;
  //   // stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   // stream.ignore(21);
  //   // stream >> m_fit.average_percent_dev;
  //   // stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  //   return true;
  // }

} // namespace corsika