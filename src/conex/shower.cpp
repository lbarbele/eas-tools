#include <utility>

#include <conex/shower.h>

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

} // namespace conex