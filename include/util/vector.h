#ifndef _utl_vector_h
#define _utl_vector_h

#include <array>
#include <ostream>
#include <iomanip>

namespace util {

  class vector : public std::array<double, 3> {
  public:
    double norm() const;
    vector& normalize(double w = 1);

    static vector cross_product(const vector& u, const vector& v);

    vector operator*(const double s) const;
    vector operator/(const double s) const;
    vector operator+(const vector& v) const;
    vector operator-(const vector& v) const;
    double operator*(const vector& v) const;

    vector& operator*=(const double s);
    vector& operator/=(const double s);
    vector& operator+=(const vector& v);
    vector& operator-=(const vector& v);
  };

  vector operator*(const double s, const vector& v);

  template<class charT>
  std::basic_ostream<charT>&
  operator<<(
    std::basic_ostream<charT>& stream,
    const vector& v
  )
  {
    using namespace std;
    int w = stream.width() == 0? 15 : stream.width();
    return stream << setw(w) << v[0] << setw(w) << v[1] << setw(w) << v[2];
  }

} // namespace util

#endif // _utl_vector_h