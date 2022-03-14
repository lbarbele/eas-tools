#include "corsika-subblock.h"

namespace corsika {

  subblock::type
  subblock::get_type()
  const
  {
    if (empty()) {
      return type::empty;
    }

    auto title = at<word_t>(0).str();

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