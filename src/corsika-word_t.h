#ifndef _corsika_word_h_
#define _corsika_word_h_

#include <cstdint>
#include <locale>
#include <string>

namespace corsika {

class word_t {
private:
  char32_t data;

public:
  explicit constexpr word_t(const char32_t c = 0) : data(c) {};

  template<class T, std::enable_if_t<sizeof(T)==sizeof(data),bool> = true>
  constexpr T& as()
  {return reinterpret_cast<T&>(data);}

  template<class T, std::enable_if_t<sizeof(T)==sizeof(data),bool> = true>
  constexpr const T& as() const
  {return reinterpret_cast<const T&>(data);}

  constexpr bool operator==(const word_t& other) const
  {return data == other.data;}
};

using word_traits = std::char_traits<word_t>;
using word_codecvt = std::codecvt<corsika::word_t, char, std::mbstate_t>;

} // namespace corsika



// specialization of the members of
// std::codecvt<corsika::word_t, char, std::mbstate_t> (aka corsika::word_codecvt)
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
) const;

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
) const;

template<>
std::codecvt_base::result
corsika::word_codecvt::do_unshift(
  state_type& state,
  extern_type* to,
  extern_type* to_end,
  extern_type*& to_next
) const;

template<> int corsika::word_codecvt::do_length(
  state_type& state,
  const extern_type* from,
  const extern_type* from_end,
  size_t max
) const;

template<>
int corsika::word_codecvt::do_encoding() const noexcept;

template<>
bool corsika::word_codecvt::do_always_noconv() const noexcept;

template<>
int corsika::word_codecvt::do_max_length() const noexcept;



// full specialization of the class std::char_traits<corsika::word_t>
// (aka corsika::word_traits)
template<>
struct std::char_traits<corsika::word_t> {
  using char_type = corsika::word_t;
  using int_type = std::char_traits<char32_t>::int_type;
  using off_type = std::char_traits<char32_t>::off_type;
  using pos_type = std::char_traits<char32_t>::pos_type;
  using state_type = std::char_traits<char32_t>::state_type;

  static constexpr void
  assign(char_type& w1, const char_type& w2)
  {w1 = w2;}

  static constexpr bool
  eq(const char_type& w1, const char_type& w2)
  {return w1 == w2;}

  static constexpr bool
  lt(const char_type& w1, const char_type& w2)
  {return w1.as<float>() < w2.as<float>();}

  static constexpr int
  compare(const char_type* s1, const char_type* s2, size_t n)
  {
    for(size_t i = 0; i < n; ++i) {
      if (lt(s1[i], s2[i])) {
        return -1;
      } else if (lt(s2[i], s1[i])) {
        return 1;
      }
    }
    return 0;
  }

  static constexpr size_t
  length(const char_type* s)
  {
    size_t i = 0;
    constexpr auto ref = char_type();
    while (!eq(s[i], ref)) {
      ++i;
    }
    return i;
  }

  static constexpr const char_type*
  find(const char_type* s, size_t n, const char_type& w)
  {
    for(size_t i = 0; i < n; i++) {
      if (eq(s[i], w)) {
        return s+i;
      }
    }
    return 0;
  }

  static char_type*
  move(char_type* s1, const char_type* s2, size_t n);

  static char_type*
  copy(char_type* s1, const char_type* s2, size_t n);

  static char_type*
  assign(char_type* s, size_t n, char_type w);

  static constexpr char_type
  to_char_type(const int_type& i)
  {return corsika::word_t(static_cast<char32_t>(i));}

  static constexpr int_type
  to_int_type(const char_type& w)
  {return int_type(w.as<char32_t>());}

  static constexpr bool
  eq_int_type(const int_type& i1, const int_type& i2)
  {return i1 == i2;}

  static constexpr int_type
  eof()
  {return char_traits<char32_t>::eof();}

  static constexpr int_type
  not_eof(const int_type& i)
  {return char_traits<char32_t>::not_eof(i);}
};


#endif // _corsika_word_h_