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
#include <TROOT.h>

#include <conex/file.h>
#include <util/math.h>
#include <models/dedx_profile.h>

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
  // let ROOT compute the chi2 in parallel
  ROOT::EnableThreadSafety();
  ROOT::EnableImplicitMT();

  // so we don't have to type the namespaces everytime
  using models::dedx_profile::get_usp_l;
  using models::dedx_profile::get_usp_r;
  using models::dedx_profile::get_gh_x0;
  using models::dedx_profile::get_gh_lambda;

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
    const auto id = cxFile.get_header().get_particle();

    // functions used to create the TF1s
    auto single_gh = [=](const double* x, const double* p) {
      const double lgecal = std::log10(p[0]) - 9;
      const double params[4] = {p[0], get_gh_x0(lgecal, id), p[1], get_gh_lambda(lgecal, id)};
      return util::math::gaisser_hillas(x, params);
    };

    auto single_usp = [=](const double* x, const double* p) {
      const double lgecal = std::log10(p[0]) - 9;
      const double params[4] = {p[0], p[1], get_usp_l(lgecal, id), get_usp_r(lgecal, id)};
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