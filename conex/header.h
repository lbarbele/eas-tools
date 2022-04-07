#ifndef _conex_header_h
#define _conex_header_h

namespace conex {

  struct header {
    int    Seed1;
    int    Particle;
    double Alpha;
    double lgEmin;
    double lgEmax;
    double zMin;
    double zMax;
    float  Version;
    float  OutputVersion;
    int    HEModel;
    int    LEModel;
    float  HiLowEgy;
    float  hadCut;
    float  emCut;
    float  hadThr;
    float  muThr;
    float  emThr;
    float  haCut;
    float  muCut;
    float  elCut;
    float  gaCut;
    
    double lambdaLgE[31];
    double lambdaProton[31];
    double lambdaPion[31];
    double lambdaHelium[31];
    double lambdaNitrogen[31];
    double lambdaIron[31];

    // conex extensions
    int    resamplingMode;
    double modThreshold;
    double f19_cx;
    double f19_meson;
    double f19;
  };

} // namespace conex

#endif // _conex_header_h