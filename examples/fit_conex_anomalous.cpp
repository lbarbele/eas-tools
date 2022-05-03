#include <iostream>
#include <algorithm>
#include <random>
#include <cmath>
#include <array>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResult.h>
#include <TString.h>
#include <TSystem.h>

#include <conex/file.h>
#include <util/math.h>

enum class particle_id {p = 100, He = 400, C = 1200};

// see the description of these functions below!
double get_L(const double lgecal, const particle_id id);
double get_R(const double lgecal, const particle_id id);
double get_X0(const double lgecal, const particle_id id);
double get_lambda(const double lgecal, const particle_id id);

TGraphErrors get_profile_cut(TGraph& g, bool doFluctuate = false);

template <int NPar>
struct FitData {
  double params[NPar];
  double errors[NPar];
  double chi2;
  unsigned int status;
};

template <int NPar>
FitData<NPar>&
operator<<(
  FitData<NPar>& fitData,
  TFitResult& fitRes
)
{
  if (fitRes.NPar() != NPar) {
    throw;
  }

  for (int i = 0; i < NPar; ++i) {
    fitData.params[i] = fitRes.Parameter(i);
    fitData.errors[i] = fitRes.Error(i);
  }

  fitData.chi2 = fitRes.MinFcnValue();
  fitData.status = fitRes.Status();

  return fitData;
}

int
main(
  int argc,
  char** argv
)
{
  // the output file
  TFile file("find_anomalous.root", "recreate");

  // the output tree and its branches
  TTree tree("data", "data");

  TString fileName;
  TGraph dedx;
  FitData<2> ghSingleFit;
  FitData<2> uspSingleFit;
  FitData<4> ghDoubleFit;
  FitData<4> uspDoubleFit;
  unsigned int ninflec = 0;
  unsigned int ishower = 0;

  tree.Branch("fileName", &fileName);
  tree.Branch("dedx", &dedx);
  tree.Branch("ninflec", &ninflec, "ninflec/i");
  tree.Branch("ishower", &ishower, "ishower/i");
  tree.Branch("ghSingleFit", &ghSingleFit, "ecal/D:xmax:ecalErr:xmaxErr:chi2:status/i");
  tree.Branch("uspSingleFit", &uspSingleFit, "ecal/D:xmax:ecalErr:xmaxErr:chi2:status/i");
  tree.Branch("ghDoubleFit", &ghDoubleFit, "ecal_1/D:xmax_1:ecal_2:xmax_2:ecalErr_1:xmaxErr_1:ecalErr_2:xmaxErr_2:chi2:status/i");
  tree.Branch("uspDoubleFit", &uspDoubleFit, "ecal_1/D:xmax_1:ecal_2:xmax_2:ecalErr_1:xmaxErr_1:ecalErr_2:xmaxErr_2:chi2:status/i");

  // loop over input files
  for (int ifile = 1; ifile < argc; ++ifile) {

    // open the conex file and check
    conex::file cxFile(argv[ifile]);
    if (!cxFile.is_open()) {
      std::cout << "bad conex file " << argv[ifile] << std::endl;
      continue;
    }

    // get file name (will go to the tree)
    fileName = gSystem->BaseName(argv[ifile]);

    // show progress
    std::cout << "\r" << ifile << "/" << argc-1 << " " << argv[ifile] << std::endl;

    // get id of the primary particle (used to get the parameters)
    const auto id = static_cast<particle_id>(cxFile.get_header().get_particle());

    // functions used to create the TF1s
    auto single_gh = [=](const double* x, const double* p) {
      const double lgecal = std::log10(p[0]) - 9;
      const double params[4] = {p[0], get_X0(lgecal, id), p[1], get_lambda(lgecal, id)};
      return util::math::gaisser_hillas(x, params);
    };

    auto single_usp = [=](const double* x, const double* p) {
      const double lgecal = std::log10(p[0]) - 9;
      const double params[4] = {p[0], p[1], get_L(lgecal, id), get_R(lgecal, id)};
      return util::math::usp_function(x, params);
    };

    auto double_gh = [&](const double*x, const double* p) {
      return single_gh(x, p) + single_gh(x, p+2);
    };

    auto double_usp = [&](const double*x, const double* p) {
      return single_usp(x, p) + single_usp(x, p+2);
    };

    // loop over showers in the current file
    for (ishower = 0; ishower < cxFile.get_n_showers(); ++ishower) {
      auto shower = cxFile.get_shower(ishower);

      // get the dEdX profile from CONEX
      dedx = shower.graph_dedx();

      // get a cut of the dEdX profile and the associated "errors"
      auto profile = get_profile_cut(dedx);

      // estimates for the calorimentric energy and Xmax
      const int imax = std::max_element(profile.GetY(), profile.GetY()+profile.GetN()) - profile.GetY();
      const double ecal = profile.Integral();
      const double xmax = profile.GetX()[imax];

      // min/max depth values
      const double min = profile.GetX()[0];
      const double max = profile.GetX()[profile.GetN()-1];

      // gaisser-hillas
      TF1 ghSingleFcn("", single_gh, min, max, 2);
      ghSingleFcn.SetParameters(ecal, xmax);
      ghSingleFit << *profile.Fit(&ghSingleFcn, "SQN");

      TF1 ghDoubleFcn("", double_gh, min, max, 4);
      ghDoubleFcn.SetParameters(0.5*ecal, xmax-200, 0.5*ecal, xmax+200);
      ghDoubleFit << *profile.Fit(&ghDoubleFcn, "SQN");

      // usp
      TF1 uspSingleFcn("", single_usp, min, max, 2);
      uspSingleFcn.SetParameters(ecal, xmax);
      uspSingleFit << *profile.Fit(&uspSingleFcn, "SQN");

      TF1 uspDoubleFcn("", double_usp, min, max, 4);
      uspDoubleFcn.SetParameters(0.5*ecal, xmax-200, 0.5*ecal, xmax+200);
      uspDoubleFit << *profile.Fit(&uspDoubleFcn, "SQN");

      // count the number of inflection points on the profile
      ninflec = util::math::count_inflection_points(profile.GetY(), profile.GetN(), 5);

      // fill the output tree
      tree.Fill();
    }
  }

  file.Write();

  return 0;
}


// create a TGraphErrors with the errors on the y coordinates
// following the procedure of arXiv:1111.0504. Optionally,
// each point can be fluctuated in the same manner as in this
// reference (set doFluctuate to true)
TGraphErrors
get_profile_cut(
  TGraph& g,
  bool doFluctuate
)
{
  const double* x = g.GetX();
  const double* y = g.GetY();
  const int nx = g.GetN();

  const double sum = std::accumulate(y, y+nx, 0.0);
  const double bound = 1e-3*sum;

  const double* first = y;
  const double* last  = y+nx-1;

  // first element larger than bound
  while (*first < bound && first <= last) {
    ++first;
  }

  // last element larger than bound
  while (*last < bound && last >= first) {
    --last;
  }

  // index of the first element larger than bound and size of the new array
  const int ifirst = first - y;
  const int npt = last + 1 - first;

  // array of "fluctuations"
  double v[npt];
  std::transform(first, first+npt, v, [=](const double& n){return 1e-2*std::sqrt(sum*n);});

  // add fluctuations to the profile
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::normal_distribution<> rand(0, 1);
  
  if (doFluctuate) {
    double yrandom[npt];
    for (int i = 0; i < npt; ++i) {
      yrandom[i] = *(first+i) + v[i]*rand(gen);
    }
    return TGraphErrors(npt, x+ifirst, yrandom, nullptr, v);
  } else {
    return TGraphErrors(npt, x+ifirst, first, nullptr, v);
  }
}

// ! parametrization of GH and USP function parameters
// these parametrizations were obtained using the tool fit_conex_profiles (see the 
// cpp file for details) using 120000 showers for each primary simulated with 
// sibyll 2.3d. The energy range of the simulations is 10^17 to 10^20 eV.
double
get_L(
  const double lgecal,
  const particle_id id
)
{
  switch (id) {
    case particle_id::p:
      return 227.183 + lgecal*(7.1665 + 0.950551*lgecal);
    case particle_id::He:
      return 229.021 + lgecal*(5.90094 + 0.523577*lgecal);
    case particle_id::C:
      return 228.763 + lgecal*(5.09137 + 0.366554*lgecal);
  }
  return 0;
}

double
get_R(
  const double lgecal,
  const particle_id id
)
{
  switch (id) {
    case particle_id::p:
      return 0.256203 - lgecal*(0.0299802 - 0.00379108*lgecal);
    case particle_id::He:
      return 0.272744 - lgecal*(0.0329877 - 0.00371185*lgecal);
    case particle_id::C:
      return 0.290858 - lgecal*(0.0362762 - 0.00336855*lgecal);
  }
  return 0;
}

double
get_X0(
  const double lgecal,
  const particle_id id
)
{
  switch (id) {
    case particle_id::p:
      return -122.225 - lgecal*(69.6578 + 4.57977*lgecal);
    case particle_id::He:
      return -112.308 - lgecal*(60.9803 + 4.84916*lgecal);
    case particle_id::C:
      return -89.7777 - lgecal*(53.2338 + 7.11282*lgecal);
  }
  return 0;
}

double
get_lambda(
  const double lgecal,
  const particle_id id
)
{
  switch (id) {
    case particle_id::p:
      return 58.7702 - lgecal*(4.95382 - 0.880473*lgecal);
    case particle_id::He:
      return 62.7939 - lgecal*(5.91256 - 0.770987*lgecal);
    case particle_id::C:
      return 66.6921 - lgecal*(6.70239 - 0.625783*lgecal);
  }
  return 0;
}
