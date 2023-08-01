#ifndef _util_math_h
#define _util_math_h

#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <type_traits>

#include <units/math.h>

namespace util::math {

  using namespace units::math;

  inline double
  gaisser_hillas(
    const double x,
    const double nmax,
    const double x0,
    const double xmax,
    const double p1,
    const double p2 = 0,
    const double p3 = 0
  )
  {
    if (x <= x0) {
      return 0;
    }

    const double l = p1 + x*(p2 + x*p3);
    const double a = (xmax - x0) / l;
    const double r = (x - x0) / (xmax - x0);
    return nmax * std::exp(a * (1 - r + std::log(r)));
  }

  inline double
  gaisser_hillas_ecal(
    const double x,
    const double ecal,
    const double x0,
    const double xmax,
    const double lambda
  )
  {
    if (x <= x0) {
      return 0;
    }
    
    const double z = (x - x0) / lambda;
    const double a = (xmax - x0) / lambda;
    return (ecal/lambda) * std::exp(-std::lgamma(a+1) - z + a*std::log(z));
  }

  inline double usp_ecal(
    const double x,
    const double ecal,
    const double xmax,
    const double l,
    const double r
  )
  {
    const double invr2 = 1.0/(r*r);
    const double z = (1.0/r + (x-xmax)/l) / r;
    return z<=0? 0 : ecal*std::pow(z,invr2)*std::exp(-z)/(l*r*std::tgamma(1+invr2));
  }
  
  template <typename T>
  unsigned int
  count_inflection_points(
    const T* f,
    const unsigned int size,
    const unsigned int step
  )
  {
    unsigned int ninflec = 0;

    T dp = 0;
    T d = f[2*step] - 2*f[step] + f[0];

    for (unsigned int i = step; i < size-2*step; i+=step) {
      dp = d;
      d = f[i+2*step] - 2*f[i+step] + f[i];

      if (dp*d <= 0) {
        ++ninflec;
      }
    }

    return ninflec;
  }

	template<class F, class T, class... Ts, class R = decltype(T{}*std::declval<F>()(T{}, Ts{}...))>
	R
	romberg_integral (
    const T a,
    const T b,
    const double eps, 
    F function,
    Ts... params
  ) {
		constexpr int max_it = 30;

		T step = b - a;
    std::array<R, max_it> v;
		v[0] = 0.5*step*(function(a, params...) + function(b, params...));

		R value_before = v[0];
		int n_steps = 1;
		
		for (int k = 1; k < max_it; k++) {
			n_steps *= 2;
			step /= 2.0;
			
      v[k] = R(0.);
			for (int i = 1; i < n_steps; i+=2)
				v[k] += function(a + i*step, params...);
			v[k] = v[k]*step + 0.5*v[k-1];
			
			for (int j = k-1; j >= 0; j--)
				v[j] = v[j+1] + (v[j+1] - v[j]) / (std::pow(4,k-j) - 1.0);
			
			if (std::fabs((value_before - v[0])/value_before) < eps)
				return v[0];
			else
				value_before = v[0];
		}
		
		throw std::runtime_error("romberg integration failed");
	}

} // namespace util::math

#endif