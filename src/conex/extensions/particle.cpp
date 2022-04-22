#include <conex/extensions/particle.h>

namespace conex::extensions {

  util::vector
  particle::get_momentum()
  const
  {
    return{Px, Py, Pz};
  }

} // conex::extensions