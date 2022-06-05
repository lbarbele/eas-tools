#ifndef _util_math_h
#define _util_math_h

namespace util::math {

  double gaisser_hillas(
    const double x,
    const double ecal,
    const double x0,
    const double xmax,
    const double lambda
  );

  double gaisser_hillas(
    const double* x,
    const double* p
  );

  double double_gaisser_hillas(
    const double x,
    const double ecal,
    const double w,
    const double x01,
    const double xmax1,
    const double lambda1,
    const double x02,
    const double xmax2,
    const double lambda2
  );

  double double_gaisser_hillas(
    const double* x,
    const double* p
  );

  double gaisser_hillas_sixpar(
    const double x,
    const double ymax,
    const double x0,
    const double xmax,
    const double p1,
    const double p2,
    const double p3
  );

  double gaisser_hillas_sixpar(
    const double* x,
    const double* p
  );

  double usp_function(
    const double x,
    const double ecal,
    const double xmax,
    const double l,
    const double r
  );
  
  double usp_function(
    const double* x,
    const double* p
  );

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

} // namespace util::math

#endif