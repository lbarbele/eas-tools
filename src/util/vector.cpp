#include <util/vector.h>

#include <cmath>

namespace util {

  // - Vector normalization

  // * compute vector norm
  constexpr
  double
  vector::norm()
  const
  {
    return std::sqrt((*this)*(*this));
  }

  // * normalize vector to given value
  constexpr 
  vector&
  vector::normalize(
    double w
  )
  {
    *this *= w/norm();
    return *this;
  }

  // - Static methods

  // * compute the cross product between two vectors
  constexpr
  vector
  vector::cross_product(
  const vector& u,
  const vector& v
  )
  {
    return vector{
      u[1]*v[2] - u[2]*v[1],
      u[2]*v[0] - u[0]*v[2],
      u[0]*v[1] - u[1]*v[0]
    };
  }

  // - Vector-scalar operations

  constexpr
  vector&
  vector::operator*=(
    const double s
  )
  {
    at(0) *= s;
    at(1) *= s;
    at(2) *= s;
    return *this;
  }

  constexpr
  vector&
  vector::operator/=(
    const double s
  )
  {
    at(0) /= s;
    at(1) /= s;
    at(2) /= s;
    return *this;
  }

  constexpr
  vector
  vector::operator*(
    const double s
  ) const
  {
    vector other = *this;
    other *= s;
    return other;
  }

  constexpr
  vector
  vector::operator/(
    const double s
  ) const
  {
    vector other = *this;
    other /= s;
    return other;
  }

  // * not a class member
  constexpr
  vector
  operator*(
    const double s,
    const vector& v
  )
  {
    return v*s;
  }


  // - Vector-vector operations

  constexpr
  vector&
  vector::operator+=(
    const vector& v
  )
  {
    at(0) += v[0];
    at(1) += v[1];
    at(2) += v[2];
    return *this;
  }

  constexpr
  vector&
  vector::operator-=(
    const vector& v
  )
  {
    at(0) -= v[0];
    at(1) -= v[1];
    at(2) -= v[2];
    return *this;
  }

  constexpr
  vector
  vector::operator+(
    const vector& v
  ) const
  {
    vector other = *this;
    other += v;
    return other;
  }

  constexpr
  vector
  vector::operator-(
    const vector& v
  ) const
  {
    vector other = *this;
    other -= v;
    return other;
  }

  constexpr
  double
  vector::operator*(
    const vector&v
  ) const
  {
    return v[0]*at(0) + v[1]*at(1) + v[2]*at(2);
  }

} 