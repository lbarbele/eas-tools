#include <corsika/longfile.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <array>

namespace corsika {

  longfile::iterator::iterator()
  : m_ishower(-1),
    m_stream(nullptr)
  {
  }

  longfile::iterator::iterator(
    const std::shared_ptr<std::ifstream>& stream
  ) :
    m_stream(stream)
  {
    if (!stream || !stream->is_open()) {
      return;
    }

    m_stream->seekg(0);

    ++(*this);
  }

  bool
  longfile::iterator::read_block()
  {
    // check if stream is ok
    if (!m_stream || !m_stream->is_open() || !m_stream->good()) {
      return false;
    }

    // clear previous content
    m_profile.clear();

    // get a reference to the stream and create a string buffer
    auto& stream = *m_stream;
    std::string buffer;

    //
    // parse the first line and skip the second one
    //
    std::getline(stream, buffer);
    if (!stream.good()) {
      // unexpected EOF
      return false;
    }

    // check first statement
    if (buffer.substr(1,28) != "LONGITUDINAL DISTRIBUTION IN") {
      return false;
    }

    // get number of longitudinal steps, shower number, step size, and depth scale
    const unsigned int nsteps = std::stoul(buffer.substr(29));
    const unsigned int ishower = std::stoul(buffer.substr(79));
    const double step_size = std::stod(buffer.substr(54));
    const std::string depth_scale = buffer.substr(36,buffer.find_first_of(' ', 36)-36);

    if (nsteps <= 0 || ishower <= 0 || step_size <= 0) {
      // bad values
      return false;
    } else {
      // resize the profile according to nsteps
      m_profile.resize(nsteps);
      // tell the profile its own step size
      m_profile.set_step_size(step_size);
      // store the number of the current shower
      m_ishower = ishower;
    }

    if (depth_scale == "VERTICAL") {
      m_profile.set_slant(false);
    } else if (depth_scale == "SLANT") {
      m_profile.set_slant(true);
    } else {
      // bad value
      return false;
    }

    // skip the second line
    std::getline(stream, buffer);
    if (!stream.good()) {
      // unexpected EOF
      return false;
    }

    //
    // read the longidutinal particle distributions
    //
    for (unsigned int iline = 0; iline < nsteps; ++iline) {
      // read all columns
      for (int iprof = 0; iprof < 10; ++iprof) {
        stream >> m_profile.get(iprof).at(iline);
      }

      // consitency check of atmospheric depth
      const double& depth = m_profile.get(0).at(iline);
      if (1e-10 < std::fabs((iline+1)*step_size/depth - 1) ) {
        std::cerr << "bad depth value when reading .long file" << std::endl;
        return false;
      }

      // check if an unexpected EOF has been reached
      if (!stream.good()) {
        return false;
      }
    }

    // ignore the rest of the last line (only the \n character)
    stream.ignore(1);

    //
    // now read the header before the energy deposit profiles
    //
    std::getline(stream, buffer);
    if (!stream.good()) {
      // unexpected EOF
      return false;
    }

    if (buffer.substr(1,30) != "LONGITUDINAL ENERGY DEPOSIT IN") {
      // bad value
      return false;
    }

    // skip the line with column names
    std::getline(stream, buffer);
    if (!stream.good()) {
      // unexpected EOF
      return false;
    }

    //
    // read the energy deposit distributions
    //
    for (unsigned int iline = 0; iline < nsteps; ++iline) {
      // read all columns
      for (int iprof = 10; iprof < 20; ++iprof) {
        stream >> m_profile.get(iprof).at(iline);
      }

      // consitency check of atmospheric depth
      const double& depth = m_profile.get(10).at(iline);
      if (1e-10 < std::fabs((iline + 0.5)*step_size/depth - 1) ) {
        std::cerr << "bad depth value when reading .long file" << std::endl;
        return false;
      }

      // check if an unexpected EOF has been reached
      if (!stream.good()) {
        return false;
      }
    }

    // ignore the rest of the last line (only the \n character)
    stream.ignore(1);

    //
    // check if there is a gaisser hillas fit and read it, if so
    //

    // store the position before reading the next line
    auto pos = stream.tellg();

    // read the first word and check if there is a fit
    stream >> buffer;

    if (buffer != "FIT") {
      stream.seekg(pos);
    } else {
      // there is a fit, skip the first line containing the formula
      std::getline(stream, buffer);

      // skip also the second line, which has no useful information
      std::getline(stream, buffer);

      // read the gaisser hillas parameters
      stream.ignore(21);

      if (!stream.good()) {
        // unexpected EOF
        return false;
      }

      std::array<double, 6> params;

      for (int i = 0; i < 6; ++i) {
        stream >> params.at(i);
      }
      std::getline(stream, buffer);

      if (!stream.good()) {
        // unexpected EOF
        return false;
      }

      util::gaisser_hillas_fit fit(params);

      // read the chi2/dof
      stream.ignore(21);

      double chi2;
      stream >> chi2;
      std::getline(stream, buffer);
      
      const int ndof = nsteps - 6;
      chi2 *= ndof;

      if (!stream.good()) {
        // unexpected EOF
        return false;
      }

      // read the average percent deviation
      stream.ignore(21);

      double percent_dev;
      stream >> percent_dev;
      std::getline(stream, buffer);
    }

    return true;
  }

  const profile&
  longfile::iterator::operator*()
  const
  {
    return m_profile;
  }

  const profile*
  longfile::iterator::operator->()
  const
  {
    return &m_profile;
  }

  longfile::iterator&
  longfile::iterator::operator++()
  {
    try {
      if (!read_block()) {
        *this = iterator();
      }
    } catch(...) {
      *this = iterator();
    }

    return *this;
  }

  longfile::iterator
  longfile::iterator::operator++(
    int
  )
  {
    iterator other = *this;
    ++(*this);
    return other;
  }

  bool
  longfile::iterator::operator==(
    const iterator& other
  ) const
  {
    return !m_stream?
      other.m_stream == nullptr : (
        m_stream == other.m_stream &&
        m_ishower == other.m_ishower
      );
  }

  bool
  longfile::iterator::operator!=(
    const iterator& other
  ) const
  {
    return !((*this) == other);
  }

} // namespace corsika