
#include <iomanip>
#include <algorithm>
#include <random>
#include <cmath>
#include <array>
#include <vector>
#include <list>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TFitResult.h>
#include <TString.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <conex/file.h>
#include <util/math.h>
#include <models/dedx_profile.h>

// ! helper functions (see definition below)
TGraphErrors getProfileCut(TGraph& g, bool doFluctuate = false);

// ! Fit classes used in the main function
// ? The purpose of these classes is to encapsulate:
// ? - the fitted function
// ? - the fit method
// ? - the resulting fit data, which will go to the output tree

// ? Virtual base class for the fitter classes
// * Implements the FitData structure, to be saved in a TTree object
// * Constructs the base TF1
// * Implements a Fit method
// * Implements a GetLeaves, returning a list of leaves of the fFitData
// * Defines the virtual operator(), the function used by the TF1 base
template<unsigned int NPar>
class VFitter : public TF1 {
private:
  struct {
    double params[NPar];
    double errors[NPar];
    double chi2;
    unsigned int status;
  } fFitData;

protected:
  std::string fDrawOpt;

public:
  VFitter() : TF1("", this, &VFitter::operator(), 0., 2000., NPar)
  {}

  virtual ~VFitter() {};

  void Draw()
  {
    TF1::Draw(fDrawOpt.c_str());
  }

  auto Data()
  {
    return reinterpret_cast<void*>(&fFitData);
  }

  std::string GetLeaves()
  {
    std::string leaves;
    for (unsigned int i = 0; i < NPar; ++i) {
      leaves += GetParName(i);
      leaves += "/D:";
    }
    for (unsigned int i = 0; i < NPar; ++i) {
      leaves += GetParName(i);
      leaves += "Err/D:";
    }
    leaves += "chi2/D:";
    leaves += "status/i";
    return leaves;
  }

  TFitResultPtr Fit(
    TGraphErrors& profile,
    const double min,
    const double max,
    const std::vector<double>& params = std::vector<double>(),
    int tries = 1
  )
  {
    TFitResultPtr result;

    // reset the function range
    SetRange(min, max);

    // if parameter estimates were given, read them
    if (params.size() == NPar) {
      SetParameters(params.data());
    }

    // try to fit until status is 0 or max tries is reached
    do {
      result = profile.Fit(this, "SQN", "", min, max);
    } while(--tries > 0 && result->Status() != 0);

    // copy the fit data to the fFitData structure
    for (unsigned int ipar = 0; ipar < NPar; ++ipar) {
      fFitData.params[ipar] = GetParameter(ipar);
      fFitData.errors[ipar] = GetParError(ipar);
    }

    // copy status and chi2 do the fFitData structure
    fFitData.chi2 = result->MinFcnValue();
    fFitData.status = !result->IsValid();

    // return the fitresult
    return result;
  }

  virtual double operator()(const double* x, const double* p) = 0;
};

// ? Mixin for adding the mass dependency on the constrained fits
// * Holds the mass ID in the protected field fMassId
// * Implements SetMass method to set fMassId
class IMassDep {
protected:
  unsigned int fMassId;
public:
  void SetMass(unsigned int id) {fMassId = id;}
};

// ? Gaisser-Hillas function, uses Ecal instead of Nmax
// ? Parameters are:
// * p[0] = Ecal
// * p[1] = X0
// * p[2] = Xmax
// * p[3] = lambda
struct GHSingleFcn : public VFitter<4> {
  GHSingleFcn(){
    SetParNames("ecal", "x0", "xmax", "lambda");
    SetLineColorAlpha(kBlue, 0.7);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) {
    return util::math::gaisser_hillas(x, p);
  }
};

// ? Double Gaisser-Hillas, also using Ecal instead of Nax
// ? Parameters are:
// * p[0] = Ecal
// * p[1] = weight of the first GH function
// * p[2] = X0_1
// * p[3] = Xmax_1
// * p[4] = lambda_1
// * p[5] = X0_2
// * p[6] = Xmax_2
// * p[7] = lambda_2
struct GHDoubleFcn : public VFitter<8> {
  GHDoubleFcn(){
    SetParNames("ecal", "w", "x0_1", "xmax_1", "lambda_1", "x0_2", "xmax_2", "lambda_2");
    SetLineColorAlpha(kBlue, 0.7);
    SetLineStyle(kDashed);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) {
    return util::math::double_gaisser_hillas(x, p);
  }
};

// ? Six-parameter Gaisser-Hillas function, same as in CONEX
// ? lambda becomes lambda(X) = p1 + X*(p2 + X*p3)
// ? Parameters are:
// * p[0] = ymax (value of the function at maximum)
// * p[1] = X0
// * p[2] = Xmax (position of the maximum)
// * p[3] = p1
// * p[4] = p2
// * p[5] = p3
struct GHSixParamFcn : public VFitter<6> {
  GHSixParamFcn(){
    SetParNames("ymax", "x0", "xmax", "p1", "p2", "p3");
    SetMarkerColorAlpha(kGreen, 0.7);
    SetMarkerStyle(kOpenCircle);
    SetMarkerSize(0.5);
    SetNpx(100);
    fDrawOpt = "same p";
  }

  double operator()(const double* x, const double* p) {
    return util::math::gaisser_hillas_sixpar(x, p);
  }
};

// ? Single Gaisser-Hillas with X0 and lambda fixed
// ? These two constrained parameters are computed in terms of
// ? the primary mass (which must be set by calling SetMass)
// ? and the value of Ecal
// ? Parameters are:
// * p[0] = Ecal
// * p[1] = Xmax
struct GHSingleConstrainedFcn : public VFitter<2>, public IMassDep {
  GHSingleConstrainedFcn(){
    SetParNames("ecal", "xmax");
    SetLineColorAlpha(kRed, 0.7);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) {
    const double& ecal = p[0];
    const double& xmax = p[1];
    const double lgecal = std::log10(ecal) - 9;
    const double x0 = models::dedx_profile::get_gh_x0(lgecal, fMassId);
    const double lambda = models::dedx_profile::get_gh_lambda(lgecal, fMassId);
    return util::math::gaisser_hillas(*x, ecal, x0, xmax, lambda);
  }
};

// ? Single Gaisser-Hillas with X0 and lambda fixed
// ? These two constrained parameters are computed in terms of
// ? the primary mass (which must be set by calling SetMass)
// ? and the value of Ecal
// ? Parameters are:
// * p[0] = Ecal
// * p[1] = weight of the first GH function
// * p[2] = Xmax_1
// * p[3] = Xmax_2
struct GHDoubleConstrainedFcn : public VFitter<4>, public IMassDep {
  GHDoubleConstrainedFcn(){
    SetParNames("ecal", "w", "xmax_1", "xmax_2");
    SetLineColorAlpha(kRed, 0.7);
    SetLineStyle(kDashed);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) {
    const double lgecal1 = std::log10(p[1]*p[0]) - 9;
    const double lgecal2 = std::log10((1-p[1])*p[0]) - 9;
    const double x01 = models::dedx_profile::get_gh_x0(lgecal1, fMassId);
    const double x02 = models::dedx_profile::get_gh_x0(lgecal2, fMassId);
    const double l1 = models::dedx_profile::get_gh_lambda(lgecal1, fMassId);
    const double l2 = models::dedx_profile::get_gh_lambda(lgecal2, fMassId);
    return util::math::double_gaisser_hillas(*x, p[0], p[1], x01, p[2], l1, x02, p[3], l2);
  }
};

// ! the output canvas
class ProfileCanvas : public TCanvas {
private:
  std::string fName;
public:
  ProfileCanvas(const std::string& name) : fName(name) {
    if (!fName.empty()) {
      TCanvas::Print((fName + "[").c_str());
    }
  }

  ~ProfileCanvas() {
    if (!fName.empty()) {
      TCanvas::Print((fName + "]").c_str());
    }
  }

  void Print() {
    if (!fName.empty()) {
      TCanvas::Print(fName.c_str());
    }
  }
};

// ! main function
// ? loops over input files, which are given through stdin
// ? fits every shower in these files
// ? builts a tree containing an entry for each fitted shower with:
// ? - the corresponding file name
// ? - the index of the shower inside the file (starting at 0)
// ? - the log10(E/eV
// ? - the dEdX profile (TGraph)
// ? - the number of inflection points on the dEdX profile
// ? - fitted data for each function
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
  // possible options are:
  // --plot <filename> 
  // --out <filename>
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
  TGraph dedx;

  GHSingleFcn ghSingleFit;
  GHDoubleFcn ghDoubleFit;
  GHSingleConstrainedFcn ghSingleConstrainedFit;
  GHDoubleConstrainedFcn ghDoubleConstrainedFit;
  GHSixParamFcn ghSixParFit;

  dataTree.Branch("ifile", &ifile, "ifile/i");
  dataTree.Branch("ishower", &ishower, "ishower/i");
  dataTree.Branch("lgE", &lgE, "lgE/D");
  dataTree.Branch("dedx", &dedx);
  dataTree.Branch("ninflec", &ninflec, "ninflec/i");

  dataTree.Branch("ghSingleFit", ghSingleFit.Data(), ghSingleFit.GetLeaves().c_str());
  dataTree.Branch("ghDoubleFit", ghDoubleFit.Data(), ghDoubleFit.GetLeaves().c_str());
  dataTree.Branch("ghSingleConstrainedFit", ghSingleConstrainedFit.Data(), ghSingleConstrainedFit.GetLeaves().c_str());
  dataTree.Branch("ghDoubleConstrainedFit", ghDoubleConstrainedFit.Data(), ghDoubleConstrainedFit.GetLeaves().c_str());
  dataTree.Branch("ghSixParFit", ghSixParFit.Data(), ghSixParFit.GetLeaves().c_str());

  // * write a tree with the names of the input files
  TTree inputFilesTree("inputFiles", "inputFiles");
  TString fileName;
  inputFilesTree.Branch("fileName", &fileName);
  while (std::cin >> fileName) {
    inputFilesTree.Fill();
  }

  // * create the output canvas (only if plotFileName is not empty)
  ProfileCanvas canvas(plotFileName);

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
      const double min = profile.GetX()[0];
      const double max = profile.GetX()[profile.GetN()-1];

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
      ghSingleFit.Fit(profile, min, max, {ecal, x0, xmax, lambda}, maxTries);
      ghDoubleFit.Fit(profile, min, max, {ecal, w, x0, xmax_1, lambda, x0, xmax_2, lambda}, maxTries);
      ghSixParFit.Fit(profile, min, max, {dedxmx, x0six, xmax, p1, p2, p3}, maxTries);
      ghSingleConstrainedFit.Fit(profile, min, max, {ecal, xmax}, maxTries);
      ghDoubleConstrainedFit.Fit(profile, min, max, {ecal, w, xmax_1, xmax_2}, maxTries);

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