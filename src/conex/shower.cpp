#include <utility>

#include <conex/shower.h>
#include <util/frame.h>

namespace conex {

  TGraph
  shower::graph_dedx()
  const
  {
    TGraph g(get_nx() - 1);
    for (int i = 0; i < get_nx() - 1; ++i) {
      const double x = 0.5*(get_depths()[i] + get_depths()[i+1]);
      g.SetPoint(i, x, get_dedx()[i]);
    }
    return g;
  }

  util::vector_d
  shower::get_axis()
  const
  {
    const double theta = get_zenith_rad();
    const double phi = get_azimuth_rad();
    return util::vector_d(
      std::sin(theta) * std::cos(phi),
      std::sin(theta) * std::sin(phi),
      std::cos(theta),
      util::frame::conex_observer
    );
  }

} // namespace conex