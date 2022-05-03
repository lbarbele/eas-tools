
namespace util::math {

  double gaisser_hillas(const double x, const double ecal, const double x0, const double xmax, const double l);
  double gaisser_hillas(const double* x, const double* p);

  double usp_function(const double x, const double ecal, const double xmax, const double l, const double r);
  double usp_function(const double* x, const double* p);

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