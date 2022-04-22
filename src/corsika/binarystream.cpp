#include <corsika/binarystream.h>
#include <corsika/binarystream-iterator.h>

namespace corsika {

  // constructors: construct the base with the equivalent constructor and
  // imbue a locale containing the word_codevt conversion facet
  // if opening a file, always set binary mode
  binarystream::binarystream()
  : std::basic_fstream<word_t>()
  {
    imbue(std::locale(std::locale(), new word_codecvt));
  }

  binarystream::binarystream(
    const char* fname,
    ios_base::openmode mode
  ) :
    std::basic_fstream<word_t>(fname, mode|ios_base::binary)
  {
    imbue(std::locale(std::locale(), new word_codecvt));
  }

  binarystream::binarystream(
    const std::string& fname,
    ios_base::openmode mode
  ) :
    std::basic_fstream<word_t>(fname, mode|ios_base::binary)
  {
    imbue(std::locale(std::locale(), new word_codecvt));
  }

  binarystream::binarystream(
    const std::filesystem::path& fname,
    ios_base::openmode mode
  ) :
    std::basic_fstream<word_t>(fname, mode|ios_base::binary)
  {
    imbue(std::locale(std::locale(), new word_codecvt));
  }

  // fstream::open: call base's open, but always adding ios_base::binary
  // to the open mode
  void
  binarystream::open(
    const char* fname,
    ios_base::openmode mode
  )
  {
    std::basic_fstream<word_t>::open(fname, mode|ios_base::binary);
  }

  void
  binarystream::open(
    const std::string& fname,
    ios_base::openmode mode
  )
  {
    std::basic_fstream<word_t>::open(fname, mode|ios_base::binary);
  }

  void
  binarystream::open(
    const std::filesystem::path& fname,
    ios_base::openmode mode
  )
  {
    std::basic_fstream<word_t>::open(fname, mode|ios_base::binary);
  }

  // begin/end iterators
  binarystream::iterator
  binarystream::begin()
  {
    return binarystream::iterator(*this);
  }

  binarystream::iterator
  binarystream::end()
  {
    return binarystream::iterator();
  }

  // insertion/extraction operators
  std::basic_istream<word_t>&
  operator>>(
    std::basic_istream<word_t>& stream,
    word_t& word
  )
  {
    stream.get(word);
    return stream;
  }

  std::basic_ostream<word_t>&
  operator<<(
    std::basic_ostream<word_t>& stream,
    word_t& word
  )
  {
    stream.put(word);
    return stream;
  }

} // namespace corsika