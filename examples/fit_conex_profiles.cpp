/*
 * fit_conex_profiles.cpp
 * last update: 02-05-2022
 * 
 * This program will take as argument a list of files produced by CONEX
 * and fit the dEdX profiles of the corresponding simulations. Fits are
 * performed both using the usual Gaisser-Hillas function and also the
 * slightly different USP profile. Note that, instead of the typical
 * "Nmax" parameter, we use the calorimetric energy in both fits.
 * 
 * After performing a first fit of both functions using all their four
 * parameters, this program parametrizes the parameters L (of the USP
 * function) and X0 (of the Gaisser-Hillas function) as a function
 * of Ecal, the fitted calorimetric energy deposited by the shower
 * (which is the first parameter of both distributions).
 * 
 * Then, the program proceeds to produce fits of both functions with only
 * three parameters free, letting the other one (either L or X0) be given
 * by the value of Ecal set by the minimizer.
 * 
 * Using the fits with three parameters, the parameters R (of the USP
 * function) and lambda (of the Gaisser-Hillas function) are also
 * parametrized as a function of Ecal, which allows for a last step
 * in which fits are produced with only Xmax and Ecal free, and the
 * other two parameters being given by the value of Ecal.
 * 
 * The intent is to select, between the USP and Gaisser-Hillas functions,
 * the one that best describe the shower profiles after the parameters
 * have been parametrized.
 * 
 * The parametrizations are all performed as a second degree polynomial
 * of Ecal, that is, if x = log10(Ecal/EeV), then:
 * 
 * f(y) = [0] + x * ([1] + [2]*x)
 *
 * Syntax
 * ------
 * 
 * $ ./fit_conex_profiles conex1.root conex2.root ...
 * 
 * Output
 * ------
 * 
 * a file called "fits.root" containing:
 * 
 * - (TProfile) parameter_L: mean value of L as a function of Ecal
 * - (TProfile) parameter_R: mean value of R as a function of Ecal
 * - (TProfile) parameter_X0: mean value of X0 as a function of Ecal
 * - (TProfile) parameter_lambda: mean value of lambda as a function of Ecal
 * - (TF1) fit_L: fit of L as a function of Ecal
 * - (TF1) fit_R: fit of R as a function of Ecal
 * - (TF1) fit_X0: fit of X0 as a function of Ecal
 * - (TF1) fit_lambda: fit of lambda as a function of Ecal
 * - (TTree) fits: contains the fitted data and all the fit results
 * - histograms directory: contains histograms of chi2 and Xmax errors
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
#include <iomanip>
#include <array>
#include <random>
#include <algorithm>
#include <map>
#include <list>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TProfile.h>
#include <TCanvas.h>

#include <conex/file.h>
#include <util/math.h>

TGraphErrors get_profile_cut(TGraph& g, bool doFluctuate = false);

void four_parameter_fits(TFile& file, int argc, char** argv);
void three_parameter_fits(TFile& file);
void two_parameter_fits(TFile& file);
void first_parametrization(TFile& file);
void second_parametrization(TFile& file);
void analyze(TFile& file);

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

  // create the output file and the fits tree
  TFile file("fits.root", "recreate");

  // I: fit CONEX profiles to both gaisser-hillas and USP function
  four_parameter_fits(file, argc, argv);
  file.Write();

  // II: parametrize L and X0 as a function of Ecal (the first fit parameter)
  first_parametrization(file);

  // III: use the parametrization to build fits with three parameters
  three_parameter_fits(file);

  // IV: parametrize R and lambda as a function of Ecal using the fits with L and X0 fixed
  second_parametrization(file);

  // V: use both parametrization to fit the profiles with two parameters, as a test
  two_parameter_fits(file);

  // VI: produce histograms of chi2 so a parametrization can be choosen
  analyze(file);

  return 0;
}



void
analyze(
  TFile& file
)
{
  std::list<std::string> fitNames = {"ghFit", "ghFit_three", "ghFit_two", "uspFit", "uspFit_three", "uspFit_two"};

  // a reader for the tree with fits
  TTreeReader reader("fits", &file);

  TTreeReaderValue<double> Xmx(reader, "Xmx");

  std::map<std::string, TTreeReaderValue<TFitResult>> fitResults;
  std::map<std::string, TH1D> chi2Histograms;
  std::map<std::string, TH1D> xmaxDevHistograms;
  std::map<std::string, TProfile> xmaxDevProfiles;
  for (const auto& fitName : fitNames) {
    fitResults.try_emplace(fitName, reader, fitName.c_str());
    chi2Histograms.try_emplace(fitName, "", fitName.c_str(), 100, 0, 100);
    xmaxDevHistograms.try_emplace(fitName, "", fitName.c_str(), 500, -250, 250);
    xmaxDevProfiles.try_emplace(fitName, "", fitName.c_str(), 100, 500, 1500);
  }

  while (reader.Next()) {
    for (auto& [fitName, fitResult] : fitResults) {
      if (fitResult->Status() == 0) {
        chi2Histograms.at(fitName).Fill(fitResult->MinFcnValue());
        const double xmaxDev = fitResult->Parameter(fitName[0]=='g'? 2 : 1) - *Xmx;
        xmaxDevHistograms.at(fitName).Fill(xmaxDev);
        xmaxDevProfiles.at(fitName).Fill(*Xmx, xmaxDev);
      }
    }
  }

  file.mkdir("histograms")->cd();

  for (auto& [fitName, chi2Hist] : chi2Histograms) {
    chi2Hist.Write(("chi2_" + fitName).c_str());
  }

  for (auto& [fitName, xmaxDevHist] : xmaxDevHistograms) {
    xmaxDevHist.Write(("xmaxDev_" + fitName).c_str());
  }

  for (auto& [fitName, xmaxDevProf] : xmaxDevProfiles) {
    xmaxDevProf.Write(("xmaxDevProfile_" + fitName).c_str());
  }
}



void
second_parametrization(
  TFile& file
)
{
  const double min = -0.9;
  const double max = 1.9;

  // printout
  std::cout << "+ parametrizing lambda and R" << std::endl;

  // a reader for the tree with fits
  TTreeReader reader("fits", &file);
  TTreeReaderValue<TFitResult> uspFit_three(reader, "uspFit_three");
  TTreeReaderValue<TFitResult> ghFit_three(reader, "ghFit_three");

  // the parameter profiles
  TProfile prof_lambda("", "X_{0}", 100, min, max);
  TProfile prof_r("", "L", 100, min, max);

  // loop over entries in the tree
  while (reader.Next()) {
    if (ghFit_three->Status() == 0) {
      const double lgecal = std::log10(ghFit_three->Parameter(0)) - 9;
      prof_lambda.Fill(lgecal, ghFit_three->Parameter(3), 1.0/ghFit_three->Error(3));
    }

    // L from the usp function
    if (uspFit_three->Status() == 0) {
      const double lgecal = std::log10(uspFit_three->Parameter(0)) - 9;
      prof_r.Fill(lgecal, uspFit_three->Parameter(3), 1.0/uspFit_three->Error(3));
    }
  }

  // fit the profiles, write the profiles, write the fits
  file.cd();

  TF1 fit_r("", "[0] + x*([1] + [2]*x)", min, max);
  auto res_r = *prof_r.Fit(&fit_r, "SQN");
  prof_r.Write("parameter_R");
  fit_r.Write("fit_R");

  TF1 fit_lambda("", "[0] + x*([1] + [2]*x)", min, max);
  auto res_lambda = *prof_lambda.Fit(&fit_lambda, "SQN");
  prof_lambda.Write("parameter_lambda");
  fit_lambda.Write("fit_lambda");

  if (res_r.Status()*res_lambda.Status() != 0) {
    std::cerr << "bad lambda R parametrization!" << std::endl;
    throw;
  }
}



void
first_parametrization(
  TFile& file
)
{
  const double maxChi2 = 10;
  const double min = -0.9;
  const double max = 1.9;

  // printout
  std::cout << "+ parametrizing X0 and L" << std::endl;

  // a reader for the tree with fits
  TTreeReader reader("fits", &file);
  TTreeReaderValue<TFitResult> uspFit(reader, "uspFit");
  TTreeReaderValue<TFitResult> ghFit(reader, "ghFit");

  // the parameter profiles
  TProfile prof_x0("", "X_{0}", 100, min, max);
  TProfile prof_l("", "L", 100, min, max);

  // loop over the fits to get the parameter values
  while (reader.Next()) {
    // x0 from the Gaisser-Hillas function
    if (ghFit->Status() == 0 && ghFit->MinFcnValue() < maxChi2 && ghFit->Parameter(1) > -2000) {
      const double lgecal = std::log10(ghFit->Parameter(0)) - 9;
      prof_x0.Fill(lgecal, ghFit->Parameter(1), 1.0/ghFit->Error(1));
    }

    // L from the usp function
    if (uspFit->Status() == 0 && uspFit->MinFcnValue() < maxChi2) {
      const double lgecal = std::log10(uspFit->Parameter(0)) - 9;
      prof_l.Fill(lgecal, uspFit->Parameter(2), 1.0/uspFit->Error(2));
    }
  }

  // fit the profiles, write the profiles, write the fits
  file.cd();

  TF1 fit_l("", "[0] + x*([1] + [2]*x)", min, max);
  auto res_l = *prof_l.Fit(&fit_l, "SQN");
  prof_l.Write("parameter_L");
  fit_l.Write("fit_L");

  TF1 fit_x0("", "[0] + x*([1] + [2]*x)", min, max);
  auto res_x0 = *prof_x0.Fit(&fit_x0, "SQN");
  prof_x0.Write("parameter_X0");
  fit_x0.Write("fit_X0");

  if (res_l.Status()*res_x0.Status() != 0) {
    std::cerr << "bad x0 L parametrization!" << std::endl;
    throw;
  }
}



void
two_parameter_fits(
  TFile& file
)
{
  // get the parametrizations
  std::unique_ptr<TF1> fit_L(file.Get<TF1>("fit_L"));
  std::unique_ptr<TF1> fit_R(file.Get<TF1>("fit_R"));
  std::unique_ptr<TF1> fit_X0(file.Get<TF1>("fit_X0"));
  std::unique_ptr<TF1> fit_lambda(file.Get<TF1>("fit_lambda"));

  // get the fits tree, to which we will add a branch
  std::unique_ptr<TTree> tree(file.Get<TTree>("fits"));

  // create the new branches
  TFitResult uspFit_two;
  TFitResult ghFit_two;

  auto uspFit_two_branch = tree->Branch("uspFit_two", &uspFit_two);
  auto ghFit_two_branch = tree->Branch("ghFit_two", &ghFit_two);

  // create a reader to parse the branches already written
  TTreeReader reader(tree.get());
  TTreeReaderValue<TGraphErrors> profile(reader, "profile");

  // printout
  std::cout << "+ reprocessing fits using the second parametrization" << std::endl;

  // loop over showers
  while(reader.Next()) {

    // printout
    const auto ientry = reader.GetCurrentEntry()+1;
    std::cout << "\rshower " << ientry << "/" << reader.GetEntries() << (ientry==reader.GetEntries()? "\n" : "") << std::flush;

    // estimates for the calorimentric energy and Xmax
    const int imax = std::max_element(profile->GetY(), profile->GetY()+profile->GetN()) - profile->GetY();
    const double ecal = profile->Integral();
    const double xmax = profile->GetX()[imax];

    // min/max depth values
    const double min = profile->GetX()[0];
    const double max = profile->GetX()[profile->GetN()-1];

    // make a gaisser-hillas fit, with the calorimetric energy as a free parameter
    auto ghFunction_two = [&](double* x, double* p){
      const double lgecal = std::log10(p[0]) - 9;
      p[1] = fit_X0->Eval(lgecal);
      p[3] = fit_lambda->Eval(lgecal);
      return util::math::gaisser_hillas(x, p);
    };
    
    TF1 gh("", ghFunction_two, min, max, 4);
    gh.SetParameter(0, ecal);
    gh.FixParameter(1, 0);
    gh.SetParameter(2, xmax);
    gh.FixParameter(3, 0);
    ghFit_two = *profile->Fit(&gh, "SQN");

    // make a USP fit, with the calorimetric energy as a free parameter
    auto uspFunction_two = [&](double* x, double* p){
      const double lgecal = std::log10(p[0]) - 9;
      p[2] = fit_L->Eval(lgecal);
      p[3] = fit_R->Eval(lgecal);
      return util::math::usp_function(x, p);
    };

    TF1 usp("", uspFunction_two, min, max, 4);
    usp.SetParameter(0, ecal);
    usp.SetParameter(1, xmax);
    usp.FixParameter(2, 0);
    usp.FixParameter(3, 0);
    uspFit_two = *profile->Fit(&usp, "SQN");

    // fill the tree
    uspFit_two_branch->Fill();
    ghFit_two_branch->Fill();

  } // loop over showers

  file.Write();
}



void
three_parameter_fits(
  TFile& file
)
{
  // get the L and X0 parametrizations
  std::unique_ptr<TF1> fit_L(file.Get<TF1>("fit_L"));
  std::unique_ptr<TF1> fit_X0(file.Get<TF1>("fit_X0"));

  // get the fits tree, to which we will add a branch
  std::unique_ptr<TTree> tree(file.Get<TTree>("fits"));

  // create the new branches
  TFitResult uspFit_three;
  TFitResult ghFit_three;

  auto uspFit_three_branch = tree->Branch("uspFit_three", &uspFit_three);
  auto ghFit_three_branch = tree->Branch("ghFit_three", &ghFit_three);

  // create a reader to parse the branches already written
  TTreeReader reader(tree.get());
  TTreeReaderValue<TGraphErrors> profile(reader, "profile");

  // printout
  std::cout << "+ reprocessing fits using the first parametrization" << std::endl;

  // looop over showers
  while(reader.Next()) {

    // printout
    const auto ientry = reader.GetCurrentEntry()+1;
    std::cout << "\rshower " << ientry << "/" << reader.GetEntries() << (ientry==reader.GetEntries()? "\n" : "") << std::flush;

    // estimates for the calorimentric energy and Xmax
    const int imax = std::max_element(profile->GetY(), profile->GetY()+profile->GetN()) - profile->GetY();
    const double ecal = profile->Integral();
    const double xmax = profile->GetX()[imax];

    // min/max depth values
    const double min = profile->GetX()[0];
    const double max = profile->GetX()[profile->GetN()-1];

    // make a gaisser-hillas fit, with the calorimetric energy as a free parameter
    auto ghFunction_three = [&](double* x, double* p){
      p[1] = fit_X0->Eval(std::log10(p[0]) - 9);
      return util::math::gaisser_hillas(x, p);
    };
    
    TF1 gh("", ghFunction_three, min, max, 4);
    gh.SetParameter(0, ecal);
    gh.FixParameter(1, 0);
    gh.SetParameter(2, xmax);
    gh.SetParameter(3, 60);
    ghFit_three = *profile->Fit(&gh, "SQN");

    // make a USP fit, with the calorimetric energy as a free parameter
    auto uspFunction_three = [&](double* x, double* p){
      p[2] = fit_L->Eval(std::log10(p[0]) - 9);
      return util::math::usp_function(x, p);
    };

    TF1 usp("", uspFunction_three, min, max, 4);
    usp.SetParameter(0, ecal);
    usp.SetParameter(1, xmax);
    usp.FixParameter(2, 0);
    usp.SetParameter(3, 0.25);
    uspFit_three = *profile->Fit(&usp, "SQN");

    // fill the tree
    uspFit_three_branch->Fill();
    ghFit_three_branch->Fill();

  } // loop over showers

  file.Write();
}



void
four_parameter_fits(
  TFile& file,
  int argc,
  char** argv
)
{
  // create the fits tree
  file.cd();
  TTree tree("fits", "fits");

  // create the tree branches and the corresponding objects
  TGraphErrors profile;
  TFitResult uspFit;
  TFitResult ghFit;
  TString fileName;

  tree.Branch("fileName", &fileName);
  tree.Branch("profile", &profile);
  tree.Branch("uspFit", &uspFit);
  tree.Branch("ghFit", &ghFit);

  // data that is simply copied from the CONEX files
  TGraph dedx;
  tree.Branch("dedx", &dedx);

  double lgE, XfirstIn, Xfirst, Xmx;
  tree.Branch("lgE", &lgE, "lgE/D");
  tree.Branch("XfirstIn", &XfirstIn, "XfirstIn/D");
  tree.Branch("Xfirst", &Xfirst, "Xfirst/D");
  tree.Branch("Xmx", &Xmx, "Xmx/D");

  // printout
  std::cout << "+ processing CONEX file and performing 4-parameter fits" << std::endl;

  // loop over the input files
  for (int ifile = 1; ifile < argc; ++ifile) {

    // open the conex file and check
    conex::file cxFile(argv[ifile]);
    if (!cxFile.is_open()) {
      std::cout << "bad conex file " << argv[ifile] << std::endl;
      continue;
    }

    // show progress
    std::cout << "\r" << ifile << "/" << argc-1 << " " << argv[ifile] << (ifile == argc-1? "\n" : "") << std::flush;

    // get the file name, which will be stored into the tree
    fileName = argv[ifile];

    // loop over showers in the current file
    for (auto& shower : cxFile) {

      // get the dEdX profile from CONEX
      dedx = shower.graph_dedx();

      // get a cut of the dEdX profile and the associated "errors"
      profile = get_profile_cut(dedx);

      // estimates for the calorimentric energy and Xmax
      const int imax = std::max_element(profile.GetY(), profile.GetY()+profile.GetN()) - profile.GetY();
      const double ecal = profile.Integral();
      const double xmax = profile.GetX()[imax];

      // min/max depth values
      const double min = profile.GetX()[0];
      const double max = profile.GetX()[profile.GetN()-1];

      // make a gaisser-hillas fit, with the calorimetric energy as a free parameter
      TF1 gh("", util::math::gaisser_hillas, min, max, 4);
      gh.SetParameter(0, ecal);
      gh.SetParameter(1, -130);
      gh.SetParameter(2, xmax);
      gh.SetParameter(3, 60);
      ghFit = *profile.Fit(&gh, "SQN");

      // make a USP fit, with the calorimetric energy as a free parameter
      TF1 usp("", util::math::usp_function, min, max, 4);
      usp.SetParameter(0, ecal);
      usp.SetParameter(1, xmax);
      usp.SetParameter(2, 225);
      usp.SetParameter(3, 0.25);
      uspFit = *profile.Fit(&usp, "SQN");

      // copy CONEX data
      lgE = shower.get_lge();
      XfirstIn = shower.get_first_interaction_inelasticty();
      Xfirst = shower.get_first_interaction_depth();
      Xmx = shower.get_xmx();

      // fill the tree
      tree.Fill();

    } // loop over showers

  } // loop over the input files

  file.Write();
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