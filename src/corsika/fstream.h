#ifndef _corsika_fstream_h
#define _corsika_fstream_h

#include <fstream>
#include <filesystem>
#include <locale>
#include <string>

#include <corsika/word_t.h>

namespace corsika {

  class fstream : public std::basic_fstream<word_t> {
  private:
  public:
    class iterator;

    fstream();
    fstream(const char* fname, ios_base::openmode mode = ios_base::in);
    fstream(const std::string& fname, ios_base::openmode mode = ios_base::in);
    fstream(const std::filesystem::path& fname, ios_base::openmode mode = ios_base::in);

    void open(const char* fname, ios_base::openmode mode = ios_base::in);
    void open(const std::string& fname, ios_base::openmode mode = ios_base::in);
    void open(const std::filesystem::path& fname, ios_base::openmode mode = ios_base::in);

    iterator begin();
    iterator end();

    template<typename T = word_t>
    T get()
    {
      word_t word;
      std::basic_fstream<word_t>::get(word);
      return word.as<T>();
    }
  };

  // word insertion/extraction operators
  std::basic_istream<word_t>&
  operator>>(std::basic_istream<word_t>& stream, word_t& word);

  std::basic_ostream<word_t>&
  operator<<(std::basic_ostream<word_t>& stream, word_t& word);

} // namespace corsika

#endif