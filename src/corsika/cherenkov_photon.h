#ifndef _corsika_cherenkov_photon_h
#define _corsika_cherenkov_photon_h

#include <corsika/particle.h>

namespace corsika{

  class cherenkov_photon : public particle {
  public:
    cherenkov_photon(const float* data, const bool has_thinning);
    cherenkov_photon(const particle& part);
    cherenkov_photon(particle&& part);

    const float& n;
    const float& x;
    const float& y;
    const float& u;
    const float& v;
    const float& t;
    const float& h;
    const float& w;
  };

} // namespace corsika

#endif // _corsika_cherenkov_photon_h