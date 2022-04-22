#include "vector.h"

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

  vector&
  vector::to_obs(
    const double phi,
    const double theta
  )
  {
    *this = to_obs(*this, phi, theta);
    return *this;
  }

  vector&
  vector::from_obs(
    const double phi,
    const double theta
  )
  {
    *this = from_obs(*this, phi, theta);
    return *this;
  }

  vector
  vector::to_obs(
    const vector& v,
    const double sinp,
    const double cosp,
    const double sint,
    const double cost
  )
  {
    vector ep = v;
    vector ep1;
    ep1[2] =  ep[0];
    ep1[1] = -ep[1] * cost  -  ep[2] * sint;
    ep1[0] =  ep[1] * sint  -  ep[2] * cost;
    ep[2]  = ep1[0];
    ep[0]  = ep1[1] * cosp + ep1[2] * sinp;
    ep[1]  = ep1[1] * sinp - ep1[2] * cosp;
    return ep;
  }

  vector
  vector::to_obs(
    const vector& v,
    const double phi,
    const double theta
  )
  {
    return to_obs(v,std::sin(phi),std::cos(phi),std::sin(theta),std::cos(theta));
  }

  vector
  vector::from_obs(
    const vector& v,
    const double sinp,
    const double cosp,
    const double sint,
    const double cost
  )
  {
    vector ep = v;
    vector ep1;
    ep1[2] =  ep[2];
    ep1[1] = -ep[0] * cosp -  ep[1] * sinp;
    ep1[0] =  ep[0] * sinp -  ep[1] * cosp;
    ep[0]  = ep1[0];
    ep[1]  = ep1[1] * cost  + ep1[2] * sint;
    ep[2]  = ep1[1] * sint  - ep1[2] * cost;
    return ep;
  }

  vector
  vector::from_obs(
    const vector& v,
    const double phi,
    const double theta
  )
  {
    return from_obs(v,std::sin(phi),std::cos(phi),std::sin(theta),std::cos(theta));
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