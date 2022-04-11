#include <iostream>

#include "particle_iterator.h"

namespace corsika {

    particle_iterator::particle_iterator()
    {}

    particle_iterator::particle_iterator(
      const fstream::iterator& it
    ) :
      m_stream_it(it)
    {
      // check if the stream is good
      if (m_stream_it == fstream::iterator()) {
        *this = particle_iterator();
        return;
      }

      // check if the stream iterator is pointing to a data block
      if (m_stream_it->get_type() != subblock::type::data) {
        *this = particle_iterator();
        return;
      }

      // we got a particle block, so start reading the first particle
      m_subblock_it = m_stream_it->begin();
      ++(*this);
    }

    const particle&
    particle_iterator::operator*()
    const
    {
      return m_particle;
    }

    const particle*
    particle_iterator::operator->()
    const
    {
      return &m_particle;
    }

    particle_iterator&
    particle_iterator::operator++()
    {
      if (m_subblock_it == m_stream_it->end()) {
        ++m_stream_it;
        m_subblock_it = m_stream_it->begin();
      }

      if (m_stream_it == fstream::iterator() || m_stream_it->empty()) {
        // this should never happen, we expect the stream iterator either to
        // point to a valid subblock. eof is never reached here
        std::cerr << "corsika::particle_iterator::operator++(): bad stream iterator" << std::endl;
        throw;
      } else if (m_stream_it->get_type() != subblock::type::data) {
        // a block that is not of data type means we have reached the end of
        // particle list for this shower
        *this = particle_iterator();
      } else {
        // here we got a good stream, so read the next particle
        m_particle = particle(&m_subblock_it->as<float>(), m_stream_it.has_thinning());
        m_subblock_it += 8;
        
        // if particle fields are all zero, there is no more particles
        for (int i = 0; i < 8; ++i) {
          if (m_particle[i] != 0) {
            return *this;
          }
        }

        *this = particle_iterator();
      }

      return *this;
    }

    particle_iterator
    particle_iterator::operator++(int)
    {
      particle_iterator other(*this);
      ++(*this);
      return other;
    }

    bool
    particle_iterator::operator==(
      const particle_iterator& other)
    const
    {
      if (m_stream_it == other.m_stream_it) {
        if (m_stream_it.get_pos() < 0) {
          return true;
        }  else if (m_subblock_it-m_stream_it->begin() == other.m_subblock_it-other.m_stream_it->begin()) {
          return true;
        }
      }
      return false;
    }

    bool
    particle_iterator::operator!=(
      const particle_iterator& other
    ) const
    {
      return !this->operator==(other);
    }

} // namespace corsika