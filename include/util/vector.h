#ifndef _utl_vector_h
#define _utl_vector_h

#include <array>
#include <ostream>
#include <iomanip>

#include <util/matrix.h>

namespace util {

  class vector : public std::array<double, 3> {
  public:
    // - Vector normalization

    // * compute vector norm 
    constexpr double norm() const;

    // * normalize vector to given value
    constexpr vector& normalize(double w = 1);

    // - Static methods

    // * compute the cross product between two vectors
    constexpr static vector cross_product(const vector& u, const vector& v);

    // - Vector-scalar operations

    // * assignment multiplication/division by scalar
    constexpr vector& operator*=(const double s);
    constexpr vector& operator/=(const double s);

    // * vector-scalar product/division
    constexpr vector operator*(const double s) const;
    constexpr vector operator/(const double s) const;

    // - Vector-vector operations

    // * assigment vector sum/subtraction
    constexpr vector& operator+=(const vector& v);
    constexpr vector& operator-=(const vector& v);

    // * vector-vector sum
    constexpr vector operator+(const vector& v) const;
    constexpr vector operator-(const vector& v) const;

    // * dot (vector-vector) product
    constexpr double operator*(const vector& v) const;

    // * unary plus/minus operators
    constexpr vector operator+() const;
    constexpr vector operator-() const;

    // - Vector-matrix operations

    template <class T>
    constexpr vector& operator*=(const square_matrix<T, 3>& m)
    {
      vector before = *this;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          at(i) += m(i,j)*before.at(j);
        }
      }
      return *this;
    }

  };

  // * scalar-vector product
  constexpr vector operator*(const double s, const vector& v);

  // - Matrix-vector product

  template <class T>
  constexpr vector operator*(
    const square_matrix<T, 3>& m,
    const vector& v
  )
  {
    vector u;
    for (size_t i = 0; i < 3; ++i) {
      u.at(i) = 0;
      for (size_t j = 0; j < 3; ++j) {
        u.at(i) += m(i,j)*v.at(j);
      }
    }
    return u;
  }

  // - Print function

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

}

#endif