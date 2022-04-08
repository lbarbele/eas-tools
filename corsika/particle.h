#ifndef _corsika_particle_h
#define _corsika_particle_h

#include <array>

namespace corsika {

  class particle {
  private:
    std::array<float, 8> m_data;
  
  public:
    particle();
    particle(const float* data);

    float&
    operator[](const unsigned int pos)
    {return m_data[pos];}

    const float&
    operator[](const unsigned int pos)
    const
    {return m_data[pos];}
  };

} // namespace corsika

#endif // _corsika_particle_h