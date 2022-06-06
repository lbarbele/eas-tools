/*
 ! fit_conex_anomalous.cpp
 ! last update: 03-06-2022
 *
 * This program takes as input (from stdin) a list of CONEX simulations
 * and tries to fit the dEdX profiles of each shower in these files.
 * The fits are produced in five different variations of the usual
 * Gaisser-Hillas function. The idea is that the comparison of these
 * different fits allows for the classification of anomalous profiles.
 * The different fit functions are:
 * 
 * a) Gaisser-Hillas function
 *    Usual Gaisser-Hillas profile, but using the calorimetric energy as
 *    a free parameter, instead of the standard "Nmax". The functional
 *    form is:
 * 
 *      dEdx(x) = ecal * z^(alpha-1) * exp(-z) / (lambda*Gamma(alpha))
 *            z = (x - x0) / lambda
 *    alpha - 1 = (xmax - x0) / lambda
 * 
 * b) Constrained Gaisser-Hillas function
 *    Same as (a), but constraining X0 and lambda to average values. Those
 *    are parametrized as a function of Ecal (see models/dedx_profile.h)
 *    The parametrizations consider also the mass dependency.
 * 
 * c) Six-parameter Gaisser-Hillas function
 *    This is the same Gaisser-Hillas function used in CONEX. In this
 *    case, the amplitude parameter is dEdX_mx, not Ecal, because there
 *    is no closed form expression to our function. The functional form
 *    is:
 * 
 *      dEdx(x) = ymax *
 *                ((x-x0)/(xmax-x0))^((xmax-x0)/lambda(x)) *
 *                exp((Xmax-X)/lambda(X))
 * 
 *    lambda(x) = p1 + x*(p2 + p3*x)
 * 
 * d) Double Gaisser-Hillas
 *    This is a weighted sum of two Gaisser-Hillas functions.
 * 
 *    dEdx(x) = w*dEdx_1(x, ecal) + (1-w)*dEdx_2(x, ecal)
 * 
 * e) Constrained double Gaisser-Hillas
 *    Same as above, but with X0 and lambda constrained to their average
 *    values, just as in case (b).
 *  
 * For details on the fitting functions, see the fitter classes in fitters.h
 !
 ! Syntax, options and input format
 * 
 * The program takes as input a list of CONEX files from stdin
 * 
 * find path/to/conex/files -name 'conex*.root | fit_conex_anomalous [options]
 * 
 * where the options are:
 * 
 * --plot file.pdf : plots fitted profiles into the specified file
 * --out file.root : sets the output file name (default is profileAnalysis.root)
 * --max-showers n : stop processing after n showers have been processed
 ! 
 ! Root output files
 *
 * The main output of this program is a root file containing two trees: one with
 * the names of the files that have been processed and, the main one, with data
 * from the fitted profiles.
 ?
 ? - "inputFiles" tree:
 *
 *   This tree contains a single branch called "fileName" of type TString. One entry
 *   is produced for each file that have been processed, even if it has been skipped
 *   due to the max-showers cutoff.
 ? 
 ? - "data" tree:
 * 
 *   The tree contains one entry for each processed shower and is structured in the
 *   following branches:
 * 
 *   . ifile (unsigned int): contains the index to get the corresponding file name
 *     from the "inputFiles" tree
 * 
 *   . ishower (unsigned int): contains the shower number (starting from 0) of the
 *     shower in the input file
 * 
 *   . lgE (double): base-10 logarithm of the primary energy in eV
 * 
 *   . minDepth (double): minimum value of atmospheric depth used in the fit. That
 *     is, points of the dedx with depth below this value are not considered in
 *     the fits
 * 
 *   . maxDepth (double): maximum value of atmospheric depth used in the fit. That
 *     is, points of the dedx with depth above this value are not considered in
 *     the fits.
 * 
 *   . ninflec (unsigned int): count of inflection points in the dEdX profile
 * 
 *   . dedx (TGraph): fitted dEdX profile, a simple copy from the CONEX file
 * 
 *   . ghSingleFit (struct): status, chi2, and parameters of function (a). contains
 *     the following leaves (all double precision fp numbers):
 *     . status         . ecal           . ecalErr
 *     . chi2           . x0             . x0Err
 *                      . xmax           . xmaxErr
 *                      . lambda         . lambdaErr
 * 
 *   . ghDoubleFit (struct): status, chi2, and parameters of function (d). contains
 *     the following leaves (all double precision fp numbers):
 *     . status         . x0_1           . x0_1Err
 *     . chi2           . xmax_1         . xmax_1Err
 *     . ecal           . lambda_1       . lambda_1Err
 *     . ecalErr        . x0_2           . x0_2Err
 *     . w              . xmax_2         . xmax_2Err
 *     . wErr           . lambda_2       . lambda_2Err
 * 
 *   . ghSingleConstrainedFit (struct): status, chi2, and parameters of function (b).
 *     contains the following leaves (all double precision fp numbers):
 *     . status         . ecal           . ecalErr
 *     . chi2           . xmax           . xmaxErr
 * 
 *   . ghDoubleConstrainedFit (struct): status, chi2, and parameters of function (e).
 *     contains the following leaves (all double precision fp numbers):
 *     . status         . xmax_1         . xmax_1Err
 *     . chi2           . xmax_2         . xmax_2Err
 *                      . ecal           . ecalErr
 *                      . w              . wErr
 * 
 *   . ghSixParFit (struct): status, chi2, and parameters of function (c). contains
 *     the following leaves (all double precision fp numbers):
 *     . status         . ymax           . ecalErr
 *     . chi2           . x0             . x0Err
 *                      . xmax           . xmaxErr
 *                      . p1             . p1Err
 *                      . p2             . p2Err
 *                      . p3             . p3Err
 * 
 */
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <string>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>

#include <conex/file.h>
#include <util/math.h>
#include <models/dedx_profile.h>

#include "fitters.h"
#include "multipagecanvas.h"

// ! helper function (see definition below)
TGraphErrors getProfileCut(TGraph& g, bool doFluctuate = false);

// ! main function
// ? loops over input files, which are given through stdin
// ? fits every shower in these files
// ? fills the output trees
int
main(
  int argc,
  char** argv
)
{
  // let ROOT compute the chi2 in parallel
  ROOT::EnableThreadSafety();
  ROOT::EnableImplicitMT();

  // disable printout messages from ROOT
  gErrorIgnoreLevel = kWarning;

  // * define the max number of tries when a fit fails
  const unsigned maxTries = 10;

  // * parse the command line
  std::string plotFileName = "";
  std::string outFileName = "profileAnalysis.root";
  unsigned int maxShowers = 0;
  int iarg = 1;

  while (iarg < argc) {
    std::string opt = argv[iarg++];

    if (iarg == argc) {
      std::cerr << "missing parameter for option " << opt << std::endl;
      return 1;
    }

    std::string arg = argv[iarg++];

    if (opt == "--plot") {
      plotFileName = arg;
    } else if (opt == "--out") {
      outFileName = arg;
    } else if (opt == "--max-showers") {
      try {
        maxShowers = std::stoul(arg);
      } catch (...) {
        std::cerr << "bad argument for --max-showers" << std::endl;
        return 1;
      }
    } else {
      std::cerr << "invalid option " << opt << std::endl;
      return 1;
    }
  }

  // * check the command line
  const bool doPlots = !plotFileName.empty();
  if (doPlots && plotFileName.substr(plotFileName.size()-4) != ".pdf") {
    std::cerr << "bad plot file name " << plotFileName << std::endl;
    return 1;
  }

  if (outFileName.substr(outFileName.size()-5) != ".root") {
    std::cerr << "bad output file name " << outFileName << std::endl;
    return 1;
  }

  // * open the output file
  TFile file(outFileName.c_str(), "recreate");

  // * create the output tree and set its branches
  TTree dataTree("data", "data");

  unsigned int ifile = 0;
  unsigned int ishower = 0;
  unsigned int ninflec = 0;
  double lgE = 0;
  double minDepth = 0;
  double maxDepth = 0;
  TGraph dedx;

  GHSingleFcn ghSingleFit;
  GHDoubleFcn ghDoubleFit;
  GHSingleConstrainedFcn ghSingleConstrainedFit;
  GHDoubleConstrainedFcn ghDoubleConstrainedFit;
  GHSixParamFcn ghSixParFit;

  dataTree.Branch("ifile", &ifile, "ifile/i");
  dataTree.Branch("ishower", &ishower, "ishower/i");
  dataTree.Branch("lgE", &lgE, "lgE/D");
  dataTree.Branch("minDepth", &minDepth, "minDepth/D");
  dataTree.Branch("maxDepth", &maxDepth, "maxDepth/D");
  dataTree.Branch("ninflec", &ninflec, "ninflec/i");
  dataTree.Branch("dedx", &dedx);

  dataTree.Branch(ghSingleFit.GetName(), ghSingleFit.Data(), ghSingleFit.GetLeaves().c_str());
  dataTree.Branch(ghDoubleFit.GetName(), ghDoubleFit.Data(), ghDoubleFit.GetLeaves().c_str());
  dataTree.Branch(ghSingleConstrainedFit.GetName(), ghSingleConstrainedFit.Data(), ghSingleConstrainedFit.GetLeaves().c_str());
  dataTree.Branch(ghDoubleConstrainedFit.GetName(), ghDoubleConstrainedFit.Data(), ghDoubleConstrainedFit.GetLeaves().c_str());
  dataTree.Branch(ghSixParFit.GetName(), ghSixParFit.Data(), ghSixParFit.GetLeaves().c_str());

  // * write a tree with the names of the input files
  TTree inputFilesTree("inputFiles", "inputFiles");
  TString fileName;
  inputFilesTree.Branch("fileName", &fileName);
  while (std::cin >> fileName) {
    inputFilesTree.Fill();
  }

  // * create the output canvas (only if plotFileName is not empty)
  MultiPageCanvas canvas(plotFileName);

  // * loop over input files
  unsigned int totalShowers = 0;
  for (ifile = 0; ifile < inputFilesTree.GetEntries(); ++ifile) {
    inputFilesTree.GetEntry(ifile);

    // try to open the conex file (check below)
    conex::file cxFile(fileName.Data());

    // print progress / check if file is open (skip if file is bad)
    std::cout << "\r" << ifile+1 << "/" << inputFilesTree.GetEntries() << " " << fileName;
    if (maxShowers != 0 && totalShowers >= maxShowers) {
      std::cout << " (skipped)" << std::endl;
      continue;
    } else if (!cxFile.is_open()) {
      std::cout << " (bad)" << std::endl;
      continue;
    } else {
      std::cout << std::endl;
    }

    // get id of the primary particle (used to get the parameters)
    const auto id = cxFile.get_header().get_particle();

    // set mass ID for the constrained fits
    ghSingleConstrainedFit.SetMass(id);
    ghDoubleConstrainedFit.SetMass(id);

    // loop over showers in the current file
    for (ishower = 0; ishower < cxFile.get_n_showers() && (maxShowers == 0 || totalShowers < maxShowers); ++ishower, ++totalShowers) {
      auto shower = cxFile.get_shower(ishower);

      // * data that will be forwarded to the output tree
      // the dEdX profile from CONEX
      dedx = shower.graph_dedx();
      // log10 of the primary energy in eV
      lgE = shower.get_lge();

      // * get a cut of the dEdX profile and the associated "errors"
      auto profile = getProfileCut(dedx);

      // * count the number of inflection points on the profile
      ninflec = util::math::count_inflection_points(profile.GetY(), profile.GetN(), 5);

      // * min/max depth values
      minDepth = profile.GetX()[0];
      maxDepth = profile.GetX()[profile.GetN()-1];

      // * parameter estimates
      // calorimetric energy: use integral of the dEdX profile
      const double ecal = dedx.Integral();
      const double lgecal = std::log10(ecal) - 9;
      // xmax: read from CONEX's six-parameter fit
      const double xmax = shower.get_xmax();
      // x0 and lambda: use our parametrization in terms of ecal
      const double x0 = models::dedx_profile::get_gh_x0(lgecal, id);
      const double lambda = models::dedx_profile::get_gh_lambda(lgecal, id);
      // parameters for the six-param. fit: read from CONEX
      const double x0six = shower.get_x0();
      const double dedxmx = shower.get_dedx_mx();
      const double p1 = shower.get_p1();
      const double p2 = shower.get_p2();
      const double p3 = shower.get_p3();
      // weight of the double Gaisser-Hillas: use 1/2
      const double w = 0.5;
      // xmax of the double Gaisser-Hillas: use +- 200 g/cm^2 around xmax
      const double xmax_1 = xmax-200;
      const double xmax_2 = xmax+200;

      // * perform the fits
      ghSingleFit.Fit(profile, minDepth, maxDepth, {ecal, x0, xmax, lambda}, maxTries);
      ghDoubleFit.Fit(profile, minDepth, maxDepth, {ecal, w, x0, xmax_1, lambda, x0, xmax_2, lambda}, maxTries);
      ghSixParFit.Fit(profile, minDepth, maxDepth, {dedxmx, x0six, xmax, p1, p2, p3}, maxTries);
      ghSingleConstrainedFit.Fit(profile, minDepth, maxDepth, {ecal, xmax}, maxTries);
      ghDoubleConstrainedFit.Fit(profile, minDepth, maxDepth, {ecal, w, xmax_1, xmax_2}, maxTries);

      // * draw
      if (doPlots) {
        TGraph graph = dedx;
        graph.SetPoint(graph.GetN(), dedx.GetX()[dedx.GetN()-1], 0);
        graph.SetPoint(graph.GetN(), 0, 0);
        graph.SetFillColorAlpha(kBlack, 0.2);
        graph.Draw("af");

        ghSingleConstrainedFit.Draw();
        ghDoubleConstrainedFit.Draw();
        ghSingleFit.Draw();
        ghDoubleFit.Draw();
        ghSixParFit.Draw();

        canvas.Print();
      }

      // * fill the output tree
      dataTree.Fill();
    }
  }

  // * write tree buffer to the output file
  file.Write();

  return 0;
}

// ? getProfileCut: create a TGraphErrors with the errors on the
// ? y coordinates following the procedure of arXiv:1111.0504.
// ? Optionally, each point can be fluctuated in the same manner
// ? as in this reference (set doFluctuate to true)
TGraphErrors
getProfileCut(
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