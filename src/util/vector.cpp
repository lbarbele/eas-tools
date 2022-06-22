#include <util/vector.h>

#include <cmath>

namespace util {

  double
  vector::norm()
  const
  {
    return std::sqrt((*this)*(*this));
  }

  vector&
  vector::normalize(
    double w
  )
  {
    *this *= w/norm();
    return *this;
  }

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

  vector
  vector::operator*(
    const double s
  ) const
  {
    return {s*this->at(0), s*this->at(1), s*this->at(2)};
  }

  vector
  vector::operator/(
    const double s
  ) const
  {
    return {this->at(0)/s, this->at(1)/s, this->at(2)/s};
  }

  vector
  vector::operator+(
    const vector& v
  ) const
  {
    return {this->at(0)+v[0], this->at(1)+v[1], this->at(2)+v[2]};
  }

  vector
  vector::operator-(
    const vector& v
  ) const
  {
    return {this->at(0)-v[0], this->at(1)-v[1], this->at(2)-v[2]};
  }

  double
  vector::operator*(
    const vector&v
  ) const
  {
    return v[0]*this->at(0) + v[1]*this->at(1) + v[2]*this->at(2);
  }

  vector&
  vector::operator*=(
    const double s
  )
  {
    this->at(0) *= s;
    this->at(1) *= s;
    this->at(2) *= s;
    return *this;
  }

  vector&
  vector::operator/=(
    const double s
  )
  {
    this->at(0) /= s;
    this->at(1) /= s;
    this->at(2) /= s;
    return *this;
  }

  vector&
  vector::operator+=(
    const vector& v
  )
  {
    this->at(0) += v[0];
    this->at(1) += v[1];
    this->at(2) += v[2];
    return *this;
  }

  vector&
  vector::operator-=(
    const vector& v
  )
  {
    this->at(0) -= v[0];
    this->at(1) -= v[1];
    this->at(2) -= v[2];
    return *this;
  }

  // not a class member:
  vector
  operator*(
    const double s,
    const vector& v
  )
  {
    return v*s;
  }
  

} // namespace util