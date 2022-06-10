#ifndef _examples_double_bump_fitters_h
#define _examples_double_bump_fitters_h

#include <array>
#include <vector>
#include <string>

#include <util/math.h>
#include <models/dedx_profile.h>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

// ! Fit classes used by fit_conex_anomalous.cpp
// ? The purpose of these classes is to encapsulate:
// ? - the fitted function
// ? - the fit method
// ? - the resulting fit data, which will go to the output tree

// ? FitData: struct to hold the fit data
// * Complains with the requirements for a POD type, so it can
// * be used to create a TBranch of a TTree
template<unsigned int NPar>
struct FitData {
  double params[NPar];
  double errors[NPar];
  double chi2;
  unsigned int status;

  void
  Fill(
    const TFitResult& fitResult
  )
  {
    if (fitResult.NPar() != NPar) {
      throw("FitData::Fill: bad number of parameters in TFitResult");
    }

    chi2 = fitResult.MinFcnValue();
    status = !fitResult.IsValid();

    for (unsigned int ipar = 0; ipar < NPar; ++ipar) {
      params[ipar] = fitResult.Parameters()[ipar];
      errors[ipar] = fitResult.Errors()[ipar];
    }
  }
};

// ? VFitter: virtual base class for the fitter classes
// * Holds a FitData structure, to be saved in a TTree object
// * Constructs the base TF1 using operator() as the corresponding function
// * Implements a custom Fit method
// * Implements a GetLeaves method, returning a list of leaves of the fFitData
// * Implements a GetColumnNames method, compatible with RDataFrame methods
// * Defines the virtual operator(), the function used by the TF1 base
template<unsigned int NPar>
class VFitter : public TF1 {
protected:
  std::string fDrawOpt;
  FitData<NPar> fFitData;

public:
  VFitter() : TF1("", this, &VFitter::operator(), 0., 2000., NPar) {}
  virtual ~VFitter() {};

  void
  Draw(std::string opt = "")
  {
    TF1::Draw((opt.empty()? fDrawOpt : opt).c_str());
  }

  FitData<NPar>*
  Data()
  {
    return &fFitData;
  }

  std::vector<std::string>
  GetColumnNames()
  const
  {
    std::vector<std::string> columns;
    for (unsigned int ipar = 0; ipar < NPar; ++ipar) {
      columns.push_back(GetName() + std::string(".") + GetParName(ipar));
    }
    for (unsigned int ipar = 0; ipar < NPar; ++ipar) {
      columns.push_back(GetName() + std::string(".") + GetParName(ipar) + std::string("Err"));
    }
    columns.push_back((GetName() + std::string(".") + "chi2"));
    columns.push_back((GetName() + std::string(".") + "status"));
    return columns;
  }

  std::string
  GetLeaves()
  const
  {
    //! warning: the output of this function must conform with the declaration
    //! of the FitData struct. Fields must be given here in the very same
    //! order as they are declared there.
    std::string leaves = "";
    for (unsigned int i = 0; i < NPar; ++i) {
      leaves += std::string(GetParName(i)) + "/D:"; 
    }
    for (unsigned int i = 0; i < NPar; ++i) {
      leaves += std::string(GetParName(i)) + "Err/D:"; 
    }
    leaves += "chi2/D:status/i";
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
    fFitData.Fill(*result);

    // return the fit result
    return result;
  }

  virtual double operator()(const double* x, const double* p) const = 0;
};

// ? IMassDep: mixin for adding the mass dependency on the constrained fits
// * Holds the mass ID in the protected field fMassId
// * Implements SetMass method to set fMassId
class IMassDep {
protected:
  unsigned int fMassId;
public:
  void SetMass(unsigned int id) {fMassId = id;}
};

// ! Implementations

// ? GHSingleFcn: Gaisser-Hillas function. uses Ecal instead of Nmax
// * p[0] = ecal
// * p[1] = x0
// * p[2] = xmax
// * p[3] = lambda
struct GHSingleFcn : public VFitter<4> {
  GHSingleFcn(){
    SetName("ghSingleFit");
    SetParNames("ecal", "x0", "xmax", "lambda");
    SetLineColorAlpha(kBlue, 0.7);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) const {
    return util::math::gaisser_hillas(x, p);
  }
};

// ? GHDoubleFcn: double Gaisser-Hillas, also using Ecal instead of Nax
// * p[0] = ecal
// * p[1] = w (weight of the first GH function)
// * p[2] = x0_1
// * p[3] = xmax_1
// * p[4] = lambda_1
// * p[5] = x0_2
// * p[6] = xmax_2
// * p[7] = lambda_2
struct GHDoubleFcn : public VFitter<8> {
  GHDoubleFcn(){
    SetName("ghDoubleFit");
    SetParNames("ecal", "w", "x0_1", "xmax_1", "lambda_1", "x0_2", "xmax_2", "lambda_2");
    SetLineColorAlpha(kBlue, 0.7);
    SetLineStyle(kDashed);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) const {
    return util::math::double_gaisser_hillas(x, p);
  }
};

// ? GHSixParamFcn: Six-parameter Gaisser-Hillas function, same as in CONEX
// ? lambda becomes lambda(X) = p1 + X*(p2 + X*p3)
// * p[0] = ymax (value of the function at maximum)
// * p[1] = x0
// * p[2] = xmax (position of the maximum)
// * p[3] = p1
// * p[4] = p2
// * p[5] = p3
struct GHSixParamFcn : public VFitter<6> {
  GHSixParamFcn(){
    SetName("ghSixParFit");
    SetParNames("ymax", "x0", "xmax", "p1", "p2", "p3");
    SetMarkerColorAlpha(kGreen, 0.7);
    SetMarkerStyle(kOpenCircle);
    SetMarkerSize(0.5);
    SetNpx(100);
    fDrawOpt = "same p";
  }

  double operator()(const double* x, const double* p) const {
    return util::math::gaisser_hillas_sixpar(x, p);
  }
};

// ? GHSingleConstrainedFcn: single Gaisser-Hillas with X0 and lambda fixed
// ? These two constrained parameters are computed in terms of
// ? the primary mass (which must be set by calling SetMass)
// ? and the value of Ecal
// ? Parameters are:
// * p[0] = ecal
// * p[1] = xmax
struct GHSingleConstrainedFcn : public VFitter<2>, public IMassDep {
  GHSingleConstrainedFcn(){
    SetName("ghSingleConstrainedFit");
    SetParNames("ecal", "xmax");
    SetLineColorAlpha(kRed, 0.7);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) const {
    const double& ecal = p[0];
    const double& xmax = p[1];
    const double lgecal = std::log10(ecal) - 9;
    const double x0 = models::dedx_profile::get_gh_x0(lgecal, fMassId);
    const double lambda = models::dedx_profile::get_gh_lambda(lgecal, fMassId);
    return util::math::gaisser_hillas(*x, ecal, x0, xmax, lambda);
  }
};

// ? GHDoubleConstrainedFcn: single Gaisser-Hillas with X0 and lambda fixed
// ? These two constrained parameters are computed in terms of
// ? the primary mass (which must be set by calling SetMass)
// ? and the value of Ecal
// ? Parameters are:
// * p[0] = ecal
// * p[1] = w (weight of the first GH function)
// * p[2] = Xmax_1
// * p[3] = Xmax_2
struct GHDoubleConstrainedFcn : public VFitter<4>, public IMassDep {
  GHDoubleConstrainedFcn(){
    SetName("ghDoubleConstrainedFit");
    SetParNames("ecal", "w", "xmax_1", "xmax_2");
    SetLineColorAlpha(kRed, 0.7);
    SetLineStyle(kDashed);
    fDrawOpt = "same l";
  }

  double operator()(const double* x, const double* p) const {
    const double lgecal1 = std::log10(p[1]*p[0]) - 9;
    const double lgecal2 = std::log10((1-p[1])*p[0]) - 9;
    const double x01 = models::dedx_profile::get_gh_x0(lgecal1, fMassId);
    const double x02 = models::dedx_profile::get_gh_x0(lgecal2, fMassId);
    const double l1 = models::dedx_profile::get_gh_lambda(lgecal1, fMassId);
    const double l2 = models::dedx_profile::get_gh_lambda(lgecal2, fMassId);
    return util::math::double_gaisser_hillas(*x, p[0], p[1], x01, p[2], l1, x02, p[3], l2);
  }
};


#endif