#include <string>

#include "subblock.h"

namespace corsika {

  std::string_view
  subblock::get_title()
  const
  {
    return at<word_t>(0).str();
  }

  subblock::type
  subblock::get_type()
  const
  {
    if (empty()) {
      return type::empty;
    }

    auto title = get_title();

    if (title == "RUNH") {
      return type::run_header;
    } else if (title == "RUNE") {
      return type::run_end;
    } else if (title == "EVTH") {
      return type::event_header;
    } else if (title == "EVTE") {
      return type::event_end;
    } else if (title == "LONG") {
      return type::longitudinal;
    } else {
      return type::data;
    }
  }

} // namespace corsika