#ifndef _corsika_subblock_h
#define _corsika_subblock_h

#include <vector>

#include "corsika-word_t.h"

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
    type get_type() const;
  };

} // namespace corsika

#endif // _corsika_subblock_h