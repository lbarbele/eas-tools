
namespace util::math {

  template <typename T>
  T
  gaisser_hillas(
    const T x,
    const T nmx,
    const T x0,
    const T xmx,
    const T l
  )
  {
    return x <= x0 ? 0 : nmx*std::pow((x-x0)/(xmx-x0),(xmx-x0)/l)*std::exp((xmx-x)/l);
  }

  template <typename T>
  T
  gaisser_hillas(const T* x, const T* p)
  {
    return gaisser_hillas(*x, p[0], p[1], p[2], p[3]);
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

} // namespace util::math