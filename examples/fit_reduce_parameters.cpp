/*
 * fit_reduce_parameters.cpp
 * last update: 30-04-2022
 * 
 * This program uses as input a file produced with fit_conex_profiles.cpp to perform
 * two things:
 * 
 * - first, it parametrizes the parameters L, R, X0, and lambda of the USP and Gaisser-
 *   Hillas fits as a function of Ecal, which is also a fit parameter of both functions
 * - then, it uses these parametrizations to fit all dEdX profiles in the input file
 *   to both functions fixing either of these parametrized parameters
 * 
 * The output file, then, will contain the following objects:
 * 
 * - two TH2D histograms with the correlation between L and R and between X0 and lambda
 * - TProfiles of L, R, X0, and lambda as a function of Ecal
 * - TF1 functions used to fit these TProfiles to the formula [0] + x*([1] + x*[2]),
 *   where x = log10(Ecal/EeV)
 * 
 * Syntax: fit_conex_profiles fits.root (fits.root was produced with fit_conex_profiles.cpp)
 * Output: a file called "reduced_fits.root"
 * 
**/

#include <iostream>
#include <memory>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TString.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

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

  // do not write histograms automatically
  TH1::AddDirectory(false);

  // chi2 < 10 for more than 99% of the simulated showers
  // when fitting with 4 parameters
  const double maxChi2 = 10;
  
  // open the input file and read the "fits" tree
  TFile fitsFile(argv[1], "read");

  TTreeReader reader("fits", &fitsFile);
  TTreeReaderValue<TGraphErrors> profile(reader, "profile");
  TTreeReaderValue<TFitResult> uspFit(reader, "uspFit");
  TTreeReaderValue<TFitResult> ghFit(reader, "ghFit");
  TTreeReaderValue<TString> fileName(reader, "fileName");

  // create an output file
  TFile outFile("reduced_fits.root", "recreate");

  // min/max values of log10(Ecal/EeV) used in the fits of the parameters
  const double minlge = -0.9;
  const double maxlge = 1.9;

  // profiles of the parameters: average value vs. fitted log10(Ecal/EeV)
  TProfile prof_l("L_profile", "L", 100, minlge, maxlge);
  TProfile prof_r("R_profile", "R", 100, minlge, maxlge);
  TProfile prof_x0("X0_profile", "X_{0}", 100, minlge, maxlge);
  TProfile prof_lambda("lambda_profile", "#lambda", 100, minlge, maxlge);

  // histograms of the correlation between the parameters
  TH2D ghCorrelation("X0_lambda_correlation", "", 100, -2000, 500, 100, 20, 120);
  TH2D uspCorrelation("L_R_correlation", "", 100, 150, 350, 100, 0.05, 0.55);

  // helper structure to hold a parameter value + its error
  struct Parameter {
    double value;
    double error;
  };

  // loop over the 4 parameter fits to fill the profiles
  for (const auto& entry : reader) {
    // gaisser hillas parameters
    if (ghFit->Status() == 0 && ghFit->MinFcnValue() < maxChi2) {
      const double ecal = ghFit->Parameter(0);
      const double lgecal = std::log10(ecal) - 9;

      Parameter x0 = {ghFit->Parameter(1), ghFit->Error(1)};
      Parameter lambda = {ghFit->Parameter(3), ghFit->Error(3)};

      if (x0.value > -2000) {
        prof_x0.Fill(lgecal, x0.value, 1.0/x0.error);
        ghCorrelation.Fill(x0.value, lambda.value, 1.0/std::hypot(x0.error, lambda.error));
      }
      prof_lambda.Fill(lgecal, lambda.value, 1.0/lambda.error);
    }

    // usp parameters
    if (uspFit->Status() == 0 && uspFit->MinFcnValue() < maxChi2) {
      const double ecal = uspFit->Parameter(0);
      const double lgecal = std::log10(ecal) - 9;

      Parameter l = {uspFit->Parameter(2), uspFit->Error(2)};
      Parameter r = {uspFit->Parameter(3), uspFit->Error(3)};

      prof_l.Fill(lgecal, l.value, 1.0/l.error);
      prof_r.Fill(lgecal, r.value, 1.0/r.error);
      uspCorrelation.Fill(l.value, r.value, 1.0/std::hypot(l.error, r.error));
    }
  }

  // parametrize the parameters as second degree pol. in log10(E/EeV)
  TF1 fit_l("L_fit", "[0] + x*([1]+[2]*x)", minlge, maxlge);
  TF1 fit_r("R_fit", "[0] + x*([1]+[2]*x)", minlge, maxlge);
  TF1 fit_x0("X0_Fit", "[0] + x*([1]+[2]*x)", minlge, maxlge);
  TF1 fit_lambda("lambda_fit", "[0] + x*([1]+[2]*x)", minlge, maxlge);

  prof_l.Fit(&fit_l, "QN");
  prof_r.Fit(&fit_r, "QN");
  prof_x0.Fit(&fit_x0, "QN");
  prof_lambda.Fit(&fit_lambda, "QN");

  // write parameter profiles and their fits to the output file
  outFile.cd();
  prof_l.Write();
  fit_l.Write();
  prof_r.Write();
  fit_r.Write();
  prof_x0.Write();
  fit_x0.Write();
  prof_lambda.Write();
  fit_lambda.Write();

  // write the correlation histograms
  ghCorrelation.Write();
  uspCorrelation.Write();

  // a tree with the reduced fits and its branches
  outFile.cd();
  TTree fits("fits", "fits");

  TFitResult uspFit_lFix;
  TFitResult uspFit_rFix;
  TFitResult ghFit_x0Fix;
  TFitResult ghFit_lambdaFix;

  fits.Branch("uspFit_lFix", &uspFit_lFix);
  fits.Branch("uspFit_rFix", &uspFit_rFix);
  fits.Branch("ghFit_x0Fix", &ghFit_x0Fix);
  fits.Branch("ghFit_lambdaFix", &ghFit_lambdaFix);

  // now fit the showers using a reduced version of both functions
  for (const auto& entry : reader) {

    // helpful printout
    std::cout
      << "\rprocessing shower " << entry+1 << " of " << reader.GetEntries()
      << (entry+1 == reader.GetEntries()? "\n" : "")
      << std::flush;

    // estimates for the calorimentric energy and Xmax
    const int imax = std::max_element(profile->GetY(), profile->GetY()+profile->GetN()) - profile->GetY();
    const double ecal = profile->Integral();
    const double xmax = profile->GetX()[imax];
    const double lgecal = std::log10(ecal) - 9;

    // min/max depth values
    const double min = profile->GetX()[0];
    const double max = profile->GetX()[profile->GetN()-1];

    // USP function with L fixed
    auto usp_lFix = [&](double* x, double *p) {
      p[2] = fit_l.Eval(std::log10(p[0]) - 9);
      return uspFunction(x, p);
    };

    TF1 fcnUsp_lFix("", usp_lFix, min, max, 4);
    fcnUsp_lFix.SetParameter(0, ecal);
    fcnUsp_lFix.SetParameter(1, xmax);
    fcnUsp_lFix.FixParameter(2, 0);
    fcnUsp_lFix.SetParameter(3, 0.25);

    // USP function with R fixed
    auto usp_rFix = [&](double* x, double *p) {
      p[3] = fit_r.Eval(std::log10(p[0]) - 9);
      return uspFunction(x, p);
    };

    TF1 fcnUsp_rFix("", usp_rFix, min, max, 4);
    fcnUsp_rFix.SetParameter(0, ecal);
    fcnUsp_rFix.SetParameter(1, xmax);
    fcnUsp_rFix.SetParameter(2, fit_l.Eval(lgecal));
    fcnUsp_rFix.FixParameter(3, 0);

    // gaisser hillas function with X0 fixed
    auto gh_x0Fix = [&](double* x, double *p) {
      p[1] = fit_x0.Eval(std::log10(p[0]) - 9);
      return ghFunction(x, p);
    };

    TF1 fcnGh_x0Fix("", gh_x0Fix, min, max, 4);
    fcnGh_x0Fix.SetParameter(0, ecal);
    fcnGh_x0Fix.FixParameter(1, 0);
    fcnGh_x0Fix.SetParameter(2, xmax);
    fcnGh_x0Fix.SetParameter(3, fit_lambda.Eval(lgecal));

    // gaisser hillas function with lambda fixed
    auto gh_lambdaFix = [&](double* x, double *p) {
      p[3] = fit_lambda.Eval(std::log10(p[0]) - 9);
      return ghFunction(x, p);
    };

    TF1 fcnGh_lambdaFix("", gh_lambdaFix, min, max, 4);
    fcnGh_lambdaFix.SetParameter(0, ecal);
    fcnGh_lambdaFix.SetParameter(1, fit_x0.Eval(lgecal));
    fcnGh_lambdaFix.SetParameter(2, xmax);
    fcnGh_lambdaFix.FixParameter(3, 0);

    // actually perform the fits
    uspFit_lFix = *profile->Fit(&fcnUsp_lFix, "SQN");
    uspFit_rFix = *profile->Fit(&fcnUsp_rFix, "SQN");
    ghFit_x0Fix = *profile->Fit(&fcnGh_x0Fix, "SQN");
    ghFit_lambdaFix = *profile->Fit(&fcnGh_lambdaFix, "SQN");

    // fill the tree
    fits.Fill();
  }

  outFile.Write();

  return 0;
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

  const double z = (invr2 + (x[0]-xmax)/(l*r));

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