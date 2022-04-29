/*
 * fit_conex_profiles.cpp
 * last update: 29-04-2022
 * 
 * This program will take as argument a list of files produced by CONEX
 * and fit the dEdX profiles of the corresponding simulations. Fits are
 * performed both using the usual Gaisser-Hillas function and also the
 * slightly different USP profile. Note that, instead of the typical
 * "Nmax" parameter, we use the calorimetric energy in both fits.
 * 
 * Syntax: fit_conex_profiles conex1.root conex2.root ...
 * Output: a file called "fits.root"
 * 
 * Output
 * ------
 * 
 * The output contains a single TTree object called "fits" with the
 * following branches:
 * 
 * profile (TGraphErrors): contains the data used in the fit (see Fit Procedure below)
 *    uspFit (TFitResult): USP fit result returned by the TGraph::Fit method
 *     ghFit (TFitResult): same as above for the Gaisser-Hillas fit
 *     fileName (TString): path to the corresponding CONEX file
 * 
 * 
 * Fit procedure
 * -------------
 * 
 * We follow the method described in arXiv:1111.0504 to compute the 
 * chi-square statistic in the fit. The chi^2 presented in this reference
 * has the advantage of being independent of the shower energy. Apart from
 * that, such statistic is simply a renormalization of the usual Poison
 * ansatz.
 * 
 * NOTE, however, we do not add fluctuations to the profiles.
 * 
 * Gaisser-Hillas fit
 * ------------------
 * 
 * Parameters:
 * - 0) Ecal -> corresponding to the total calorimetric energy deposit
 * - 1) X0 -> pseudo first interaction slant depth position
 * - 2) Xmax -> depth of shower maximum
 * - 3) lambda -> scale parameter
 * 
 * Formula:
 *   dEdX(X) = Ecal * z^(alpha-1) * exp(-z) / (lambda*Gamma(alpha))
 *   where
 *   - z = (X - X0) / lambda
 *   - alpha-1 = (Xmax - X0) / lambda
 * 
 * USP fit
 * -------
 * 
 * Parameters:
 * - 0) Ecal -> corresponding to the total calorimetric energy deposit
 * - 1) Xmax -> depth of shower maximum
 * - 2) L -> shape parameter
 * - 3) R -> shape parameter
 * 
 * Formula:
 *   dEdX(X) = Ecal * z^(1/R^2) * exp(-z) / (L*R*Gamma(1+1/R^2))
 *   where z = (1/R + (X - Xmax)/L) / R
 * 
**/

#include <iostream>
#include <array>
#include <random>
#include <algorithm>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>

#include <conex/file.h>

TGraph get_dedx_profile(const conex::shower& sh);
TGraphErrors get_profile_cut(TGraph& g, bool doFluctuate = false);

double ghFunction(const double* x, const double* p);
double uspFunction(const double* x, const double* p);

int
main(
  int argc,
  char** argv
)
{
  // let ROOT compute the chi2 in parallel
  ROOT::EnableThreadSafety();
  ROOT::EnableImplicitMT();
  
  // open the output file and create the "fits" tree
  TFile file("fits.root", "recreate");
  TTree tree("fits", "fits");

  // create the tree branches and the corresponding objects
  TGraphErrors profile;
  TFitResult uspFit;
  TFitResult ghFit;
  TString fileName;

  tree.Branch("profile", &profile);
  tree.Branch("uspFit", &uspFit);
  tree.Branch("ghFit", &ghFit);
  tree.Branch("fileName", &fileName);

  // loop over the input files
  for (int ifile = 1; ifile < argc; ++ifile) {

    // open the conex file and check
    conex::file cxFile(argv[ifile]);
    if (!cxFile.is_open()) {
      std::cout << "bad conex file " << argv[ifile] << std::endl;
      continue;
    }

    // show progress
    std::cout << ifile << "/" << argc-1 << " " << argv[ifile] << std::endl;

    // get the file name, which will be stored into the tree
    fileName = argv[ifile];

    // loop over showers in the current file
    for (auto& shower : cxFile) {

      // get the dEdX profile from CONEX
      auto dedxProfile = get_dedx_profile(shower);

      // get a cut of the dEdX profile and the associated "errors"
      profile = get_profile_cut(dedxProfile);

      // estimates for the calorimentric energy and Xmax
      const int imax = std::max_element(profile.GetY(), profile.GetY()+profile.GetN()) - profile.GetY();
      const double ecal = profile.Integral();
      const double xmax = profile.GetX()[imax];

      // min/max depth values
      const double min = profile.GetX()[0];
      const double max = profile.GetX()[profile.GetN()-1];

      // make a gaisser-hillas fit, with the calorimetric energy as a free parameter
      TF1 gh("", ghFunction, min, max, 4);
      gh.SetParameter(0, ecal);
      gh.SetParameter(1, -130);
      gh.SetParameter(2, xmax);
      gh.SetParameter(3, 60);

      ghFit = *profile.Fit(&gh, "SQN");

      // make a USP fit, with the calorimetric energy as a free parameter
      TF1 usp("", uspFunction, min, max, 4);
      usp.SetParameter(0, ecal);
      usp.SetParameter(1, xmax);
      usp.SetParameter(2, 225);
      usp.SetParameter(3, 0.25);

      uspFit = *profile.Fit(&usp, "SQN");

      // fill the tree
      tree.Fill();

    } // loop over showers

  } // loop over the input files

  file.Write();

  return 0;
}

// build a TGraph with the dEdX profile for the shower
TGraph
get_dedx_profile(
  const conex::shower& sh
)
{
  TGraph g(sh.get_nx() - 1);
  for (int i = 0; i < sh.get_nx() - 1; ++i) {
    const double x = 0.5*(sh.get_depths()[i] + sh.get_depths()[i+1]);
    g.SetPoint(i, x, sh.get_dedx()[i]);
  }
  return g;
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

// universal shower profile function
double
uspFunction(
  const double *x,
  const double *p
)
{
  const double& ecal = p[0];
  const double& xmax = p[1];
  const double& l = p[2];
  const double& r = p[3];
  const double invr2 = 1.0/(r*r);

  const double z = (1.0/r + (x[0]-xmax)/l) / r;

  return z<=0? 0 : ecal*std::pow(z,invr2)*std::exp(-z)/(l*r*std::tgamma(1+invr2));
}

// gaisser hillas function
double
ghFunction(
  const double* x,
  const double* p
)
{
  const double& ecal = p[0];
  const double& x0 = p[1];
  const double& xmax = p[2];
  const double& l = p[3];

  const double z = (*x - x0)/l;

  if (z <= 0) {
    return 0;
  }

  const double am1 = (xmax - x0)/l;

  return ecal * std::pow(z, am1) * std::exp(-z) / (l * std::tgamma(1+am1));
}