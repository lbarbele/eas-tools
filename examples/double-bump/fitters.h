#ifndef _examples_double_bump_fitters_h
#define _examples_double_bump_fitters_h

#include <util/math.h>
#include <models/dedx_profile.h>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>

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


#endif