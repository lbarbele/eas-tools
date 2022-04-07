#include "corsika-word_t.h"

#include <cstdint>
#include <locale>
#include <string>

//
// impose the size of word_t to be of 4 bytes
//
static_assert(sizeof(corsika::word_t) == 4);

//
// definition of non-constexpr methods of std::char_traits<corsika::word_t>
//
std::char_traits<corsika::word_t>::char_type*
std::char_traits<corsika::word_t>::move(
  char_type* s1,
  const char_type* s2,
  size_t n
)
{
  return reinterpret_cast<corsika::word_t*>(
    std::char_traits<char32_t>::move(
      reinterpret_cast<char32_t*>(s1),
      reinterpret_cast<const char32_t*>(s2),
      n
    )
  );
}

std::char_traits<corsika::word_t>::char_type*
std::char_traits<corsika::word_t>::copy(
  char_type* s1,
  const char_type* s2,
  size_t n
)
{
  return reinterpret_cast<corsika::word_t*>(
    std::char_traits<char32_t>::copy(
      reinterpret_cast<char32_t*>(s1),
      reinterpret_cast<const char32_t*>(s2),
      n
    )
  );
}

std::char_traits<corsika::word_t>::char_type*
std::char_traits<corsika::word_t>::assign(
  char_type* s,
  size_t n,
  char_type w
)
{
  for(size_t i = 0; i < n; ++i) {
    s[i] = w;
  }
  return s;
}

//
// specialization of (the methods of) std::codecvt<corsika::word_t, char, std::mbstate_t>
//
template<>
std::codecvt_base::result
corsika::word_codecvt::do_in(
  state_type& state,
  const extern_type* from,
  const extern_type* from_end,
  const extern_type*& from_next,
  intern_type* to,
  intern_type* to_end,
  intern_type*& to_next
) const
{
  long char_size = sizeof(corsika::word_t);
  auto first = reinterpret_cast<const intern_type*>(from);
  auto last = first + std::min((from_end-from)/char_size, to_end-to);
  to_next = std::copy(first, last, to);
  from_next = reinterpret_cast<const extern_type*>(last);
  return (from_next == from_end) ? ok : partial;
}

template<>
std::codecvt_base::result
corsika::word_codecvt::do_out(
  state_type& state,
  const intern_type* from,
  const intern_type* from_end,
  const intern_type*& from_next,
  extern_type* to,
  extern_type* to_end,
  extern_type*& to_next
) const
{
  long char_size = sizeof(corsika::word_t);
  auto first = reinterpret_cast<const extern_type*>(from);
  auto last = first + char_size*std::min(from_end-from, (to_end-to)/char_size);
  to_next = std::copy(first, last, to);
  from_next = reinterpret_cast<const intern_type*>(last);
  return (from_next == from_end) ? ok : partial;
}

template<>
std::codecvt_base::result
corsika::word_codecvt::do_unshift(
  state_type& state,
  extern_type* to,
  extern_type* to_end,
  extern_type*& to_next
) const
{
  to_next = to;
  return result::ok;
}

template<>
int
corsika::word_codecvt::do_encoding()
const
noexcept
{return sizeof(corsika::word_t);}

template<>
bool
corsika::word_codecvt::do_always_noconv()
const
noexcept
{return false;}

template<>
int
corsika::word_codecvt::do_length(
  state_type& state,
  const extern_type* from,
  const extern_type* from_end,
  size_t max
) const
{return sizeof(corsika::word_t)*min(max, size_t(from_end-from)/sizeof(corsika::word_t));}

template<>
int
corsika::word_codecvt::do_max_length()
const
noexcept
{return sizeof(corsika::word_t);}