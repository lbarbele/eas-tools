#ifndef _corsika_subblock_h
#define _corsika_subblock_h

#include <vector>
#include <string>

#include "word_t.h"

namespace corsika {

  class subblock : public std::vector<word_t> {
  public:
    enum class type {
      run_header,
      run_end,
      event_header,
      event_end,
      longitudinal,
      data,
      empty
    };

  private:
  public:
    std::string_view get_title() const;
    type get_type() const;

    float& operator[](size_type pos)
    {return std::vector<word_t>::at(pos).as<float>();}

    const float& operator[](size_type pos) const
    {return std::vector<word_t>::at(pos).as<float>();}

    template<typename T = float>
    T& at(size_type pos)
    {return std::vector<word_t>::at(pos).as<T>();}

    template<typename T = float>
    const T& at(size_type pos) const
    {return std::vector<word_t>::at(pos).as<T>();}
  };

} // namespace corsika

#endif // _corsika_subblock_h