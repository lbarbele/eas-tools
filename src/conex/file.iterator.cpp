#include <conex/file.h>

namespace conex {

  file::iterator::iterator()
  : m_pos(-1),
    m_file(nullptr)
  {}

  file::iterator::iterator(
    file& f
  )
  : m_pos(0),
    m_file(&f)
  {
    f.get_shower(0);
  }

  const shower&
  file::iterator::operator*()
  {
    return m_file->get_shower(m_pos);
  }

  const shower*
  file::iterator::operator->()
  {
    return &m_file->get_shower(m_pos);
  }

  file::iterator&
  file::iterator::operator++()
  {
    if (++m_pos >= m_file->get_n_showers()) {
      *this = iterator();
    }
    return *this;
  }

  file::iterator
  file::iterator::operator++(
    int
  )
  {
    auto other = *this;
    ++(*this);
    return other;
  }

  bool
  file::iterator::operator==(
    const iterator& other
  ) const
  {
    return m_file == other.m_file && m_pos == other.m_pos;
  }

  bool
  file::iterator::operator!=(
    const iterator& other
  ) const
  {
    return !(*this == other);
  }

} // namespace conex