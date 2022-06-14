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
#include <iomanip>
#include <numeric>
#include <cmath>
#include <string>
#include <any>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <TString.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <TGraphErrors.h>

#include <Math/SpecFunc.h>

#include <tclap/CmdLine.h>

#include <conex/file.h>
#include <util/math.h>

// #include "fitters.h"
#include "multipagecanvas.h"

// ! fitters
class GHSixParFcn : public TF1 {
public:
  GHSixParFcn() : TF1("ghSixParFit", this, 0, 1, 6) {}

  struct {
    double chi2;
    double nmax;
    double x0;
    double xmax;
    double p1;
    double p2;
    double p3;
  } params;

  auto MakeBranch(TTree& tree) {
    return tree.Branch(GetName(), &params, "chi2/D:nmax:x0:xmax:p1:p2:p3");
  }

  auto Fit(
    TGraphErrors& profile, 
    const double min, 
    const double max,
    const conex::shower& shower,
    unsigned int maxTries = 10
  )
  {
    // * adjust the profile range 
    SetRange(min, max);

    // * set initial parameter values
    SetParameters(
      shower.get_dedx_mx(),
      shower.get_x0(),
      shower.get_xmax(),
      shower.get_p1(),
      shower.get_p2(),
      shower.get_p3()
    );

    for (int i = 0; i < GetNpar(); ++i) {
      SetParError(i, 0.01*GetParameter(i));
    }

    // * perform the fit
    TFitResultPtr result;
    int ntries = maxTries;
    do {
      result = profile.Fit(this, "SQN", "", min, max);
    } while (--ntries > 0 && result->Status() != 0);

    // * copy the parameters
    params.chi2 = result->MinFcnValue() * (result->Status() == 0? 1 : -1);
    params.nmax = GetParameter(0);
    params.x0 = GetParameter(1);
    params.xmax = GetParameter(2);
    params.p1 = GetParameter(3);
    params.p2 = GetParameter(4);
    params.p3 = GetParameter(5);

    return TFitResultPtr();
  }

  double operator()(const double* x, const double* p) const {
    return util::math::gaisser_hillas(x[0], p[0], p[1], p[2], p[3], p[4], p[5]);
  }
};

class GHSingleFcn : public TF1 {
private:
  double fMinDepth = 0;
  double fMaxDepth = 1;

  double GetEcal(const double* p) const {
    const double zm = (p[2] - p[1]) / p[3];
    const double za = (fMinDepth - p[1]) / p[3];
    const double zb = (fMaxDepth - p[1]) / p[3];
    return p[0] / (ROOT::Math::inc_gamma(zm+1, zb) - ROOT::Math::inc_gamma(zm+1, za));
  }

public:
  GHSingleFcn() : TF1("ghSingleFit", this, 0, 1, 4) {}

  struct {
    double chi2;
    double ecal;
    double x0;
    double xmax;
    double lambda;
  } params;

  auto MakeBranch(TTree& tree) {
    return tree.Branch(GetName(), &params, "chi2/D:ecal:x0:xmax:lambda");
  }

  auto Fit(
    TGraphErrors& profile, 
    const double min, 
    const double max,
    const double integral,
    const double xmax,
    unsigned int maxTries = 10
  ) {
    // * adjust the profile range 
    SetRange(min, max);
    fMinDepth = min;
    fMaxDepth = max;

    // * set guesses for x0, xmax, and lambda
    SetParameter(0, integral);
    SetParameter(1, -100);
    SetParameter(2, xmax);
    SetParameter(3, 60);

    // * perform the fit
    TFitResultPtr result;
    int ntries = maxTries;
    do {
      result = profile.Fit(this, "SQN", "", min, max);
    } while (--ntries > 0 && result->Status() != 0);

    // * copy the parameters
    params.chi2 = result->MinFcnValue() * (result->Status() == 0? 1 : -1);
    params.ecal = GetEcal(GetParameters());
    params.x0 = GetParameter(1);
    params.xmax = GetParameter(2);
    params.lambda = GetParameter(3);

    return result;
  }

  double operator()(const double* x, const double* p) const {
    const double ecal = GetEcal(p);
    return util::math::gaisser_hillas_ecal(*x, ecal, p[1], p[2], p[3]);
  }
};

class GHDoubleFcn : public TF1 {
private:
  double fMinDepth = 0;
  double fMaxDepth = 1;

  double GetEcal(const double *p) const{
    double ecal = 0;

    const double zm1 = (p[3] - p[2]) / p[4];
    const double za1 = (fMinDepth - p[2]) / p[4];
    const double zb1 = (fMaxDepth - p[2]) / p[4];
    ecal += p[1] * (ROOT::Math::inc_gamma(zm1+1, zb1) - ROOT::Math::inc_gamma(zm1+1, za1));

    const double zm2 = (p[6] - p[5]) / p[7];
    const double za2 = (p[8] - p[5]) / p[7];
    const double zb2 = (p[9] - p[5]) / p[7];
    ecal += (1-p[1]) * (ROOT::Math::inc_gamma(zm2+1, zb2) - ROOT::Math::inc_gamma(zm2+1, za2));

    ecal = p[0] / ecal;

    return ecal;
  }

public:
  GHDoubleFcn() : TF1("ghDoubleFit", this, 0, 1, 8) {}

  struct {
    double chi2;
    double ecal;
    double w;
    double x0_1;
    double xmax_1;
    double lambda_1;
    double x0_2;
    double xmax_2;
    double lambda_2;
  } params;

  auto MakeBranch(TTree& tree) {
    return tree.Branch(GetName(), &params, "chi2/D:ecal:w:x0_1:xmax_1:lambda_1:x0_2:xmax_2:lambda_2");
  }

  auto Fit(
    TGraphErrors& profile, 
    const double min, 
    const double max,
    const double integral,
    unsigned int maxTries = 10
  ) {
    // * adjust the profile range 
    SetRange(min, max);
    fMinDepth = min;
    fMaxDepth = max;

    // * set parameter guesses
    SetParameter(0, integral);
    SetParameter(1, 0.5);
    SetParameter(2, -120);
    SetParameter(3, 0.4*(max-min));
    SetParameter(4, 60);
    SetParameter(5, -120);
    SetParameter(6, 0.6*(max-min));
    SetParameter(7, 60);

    // * set parameter ranges
    SetParLimits(0, 0.5*integral, 1.5*integral);
    SetParLimits(1, 0, 1);
    SetParLimits(2, -5000, min);
    SetParLimits(3, min, max);
    SetParLimits(4, 1, 500);
    SetParLimits(5, -5000, min);
    SetParLimits(6, min, max);
    SetParLimits(7, 1, 500);

    // * perform the fit
    TFitResultPtr result;
    int ntries = maxTries;
    do {
      result = profile.Fit(this, "SQN", "", min, max);
    } while (--ntries > 0 && result->Status() != 0);

    // * copy the parameters
    params.chi2 = result->MinFcnValue() * (result->Status() == 0? 1 : -1);
    params.ecal = GetEcal(GetParameters());
    params.w = GetParameter(1);
    params.x0_1 = GetParameter(2);
    params.xmax_1 = GetParameter(3);
    params.lambda_1 = GetParameter(4);
    params.x0_2 = GetParameter(5);
    params.xmax_2 = GetParameter(6);
    params.lambda_2 = GetParameter(7);

    return result;
  }

  double operator()(const double* x, const double* p) const {
    const double ecal = GetEcal(p);
    return p[1] * util::math::gaisser_hillas_ecal(*x, ecal, p[2], p[3], p[4])
      + (1 - p[1]) * util::math::gaisser_hillas_ecal(*x, ecal, p[5], p[6], p[7]);
  }
};

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
  // * open the CONEX files piped to stdin and check
  conex::file cxFiles(std::cin);
  if (!cxFiles.is_open() || cxFiles.get_n_files() == 0) {
    std::cerr << "unable to open conex files" << std::endl;
    return 1;
  }

  // * define parameters and parse the command line
  TCLAP::CmdLine cmdLine("fit_conex_anomalous");

  TCLAP::ValueArg<unsigned int> maxShowers("m", "max-showers",
    "Maximum number of showers to be analyzed",
    false, -1, "nshowers", cmdLine);
  TCLAP::ValueArg<std::string> plotFile("p", "plot",
    "Create a pdf file at the given path with plots of the fitted profiles",
    false, "", "path", cmdLine);
  TCLAP::ValueArg<std::string> outFileName("o", "out",
    "Path to the output .root file (may be overwritten)",
    false, "profileAnalysis.root", "path", cmdLine);
  TCLAP::SwitchArg disableMT("s", "single-thread",
    "Disable ROOT's multithreading capabilities",
    cmdLine, false);

  cmdLine.parse(argc, argv);

  // * configure ROOT
  if (!disableMT) {
    ROOT::EnableThreadSafety();
    ROOT::EnableImplicitMT();
  }

  // * open the output file and check
  TFile file(outFileName.getValue().c_str(), "recreate");
  if (!file.IsOpen()) {
    std::cerr << "Failed to open the output file" << std::endl;
    return 1;
  }

  gErrorIgnoreLevel = kWarning;

  // * create the output tree and define its basic branches
  TTree dataTree("data", "data");

  TString fileName;
  dataTree.Branch("fileName", &fileName);

  unsigned int ishower = 0;
  dataTree.Branch("ishower", &ishower, "ishower/i");

  unsigned int ninflec = 0;
  dataTree.Branch("ninflec", &ninflec, "ninflec/i");
  
  double lgE = 0;
  dataTree.Branch("lgE", &lgE, "lgE/D");

  double minDepth = 0;
  dataTree.Branch("minDepth", &minDepth, "minDepth/D");

  double maxDepth = 0;
  dataTree.Branch("maxDepth", &maxDepth, "maxDepth/D");

  TGraph dedx;
  dataTree.Branch("dedx", &dedx);

  // * create the fit functions and connect them to the tree
  GHSingleFcn ghSingleFcn;
  ghSingleFcn.MakeBranch(dataTree);

  GHDoubleFcn ghDoubleFcn;
  ghDoubleFcn.MakeBranch(dataTree);

  GHSixParFcn ghSixParFcn;
  ghSixParFcn.MakeBranch(dataTree);

  // * the multiple-page canvas that will be used for plots
  MultiPageCanvas canvas(plotFile);

  // * loop over showers in the conex files
  unsigned int totalShowers = 0;
  for (const auto& shower : cxFiles) {

    // * check if file name has changed and print, if so
    auto currentFileName = gSystem->BaseName(cxFiles.get_current_file_name().c_str());
    if (fileName != currentFileName) {
      fileName = currentFileName;
      ishower = 0;
      std::cout
        << std::setw(5) << cxFiles.get_current_file_number()+1
        << "/" << cxFiles.get_n_files() << " "
        << fileName
        << std::endl;
    }

    // * increment the counters and check for maxShowers
    ++ishower;
    
    if (++totalShowers > maxShowers) {
      break;
    }

    // * the dEdX profile that will be fitted
    dedx = shower.graph_dedx();

    // * log of the primary energy in eV
    lgE = shower.get_lge();

    // * closed interval of the dEdX profile that will be fitted
    const double accum = std::accumulate(dedx.GetY(), dedx.GetY()+dedx.GetN(), 0.0);

    auto itFirst = dedx.GetY();
    while (*itFirst < 0.001*accum) {
      ++itFirst;
    }

    auto itLast = dedx.GetY() + dedx.GetN() - 1;
    while (*itLast < 0.001*accum) {
      --itLast;
    }

    unsigned iFirst = itFirst - dedx.GetY();
    unsigned iLast = itLast - dedx.GetY();

    minDepth = dedx.GetX()[iFirst];
    maxDepth = dedx.GetX()[iLast];

    // * number of inflection points in the dEdX profile
    ninflec = util::math::count_inflection_points(itFirst, 1+itLast-itFirst, 5);

    // * energy deposited in the fitted range
    double edep = 0;
    for (unsigned int i = iFirst; i < iLast; ++i) {
      edep += 0.5 * (dedx.GetX()[i+1] - dedx.GetX()[i]) * (dedx.GetY()[i+1] + dedx.GetY()[i]);
    }

    // * define weights for the chi2 fit
    TGraphErrors prof(dedx.GetN(), dedx.GetX(), dedx.GetY(), nullptr, nullptr);
    for (int i = 0; i < dedx.GetN(); ++i) {
      prof.SetPointError(i, 0, std::sqrt(1e-2 * dedx.GetY()[i] * accum));
    }

    // * perform the fits
    ghSingleFcn.Fit(prof, minDepth, maxDepth, edep, shower.get_xmx());
    ghDoubleFcn.Fit(prof, minDepth, maxDepth, edep);
    ghSixParFcn.Fit(prof, minDepth, maxDepth, shower);

    // * fill the data tree
    dataTree.Fill();

    // * draw
    if (plotFile.getValue().size() > 0) {

      dedx.SetPoint(dedx.GetN(), dedx.GetX()[dedx.GetN()-1], 0);
      dedx.SetPoint(dedx.GetN(), 0, 0);
      dedx.SetFillColorAlpha(kBlue, 0.6);
      dedx.Draw("af");

      ghDoubleFcn.SetLineColorAlpha(kRed, 0.8);
      ghDoubleFcn.SetLineWidth(3);
      ghDoubleFcn.Draw("same l");

      ghSixParFcn.SetLineColorAlpha(kGreen, 0.9);
      ghSixParFcn.SetLineStyle(kDashed);
      ghSixParFcn.Draw("same l");

      canvas.Print();
    }
  }

  // * write the tree's buffer to the output file
  file.Write();

  return 0;
}