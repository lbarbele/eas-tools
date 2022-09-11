#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>
#include <list>
#include <functional>

#include <TGraph.h>
#include <TF1.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>

#include <util/math.h>

// #include "fitters.h"
#include "multipagecanvas.h"

struct TH1DLogModel : public ROOT::RDF::TH1DModel {
  TH1DLogModel() = delete;
  TH1DLogModel(
    const char* name,
    const char* title,
    const unsigned int nbins,
    const double min,
    const double max
  ) : ROOT::RDF::TH1DModel(name, title, nbins, std::vector<double>(nbins+1, min).data())
  {
    for (unsigned int ibin = 0; ibin <= nbins; ++ibin) {
      fBinXEdges.at(ibin) *= std::pow(max/min, double(ibin)/nbins);
    }
  }
};

struct TH2DLogModel : public ROOT::RDF::TH2DModel {
  TH2DLogModel() = delete;
  TH2DLogModel(
    const char* name,
    const char* title,
    const unsigned int nbinsx,
    const double xmin,
    const double xmax,
    const unsigned int nbinsy,
    const double ymin,
    const double ymax,
    const bool logx = true,
    const bool logy = true
  ) : ROOT::RDF::TH2DModel(name, title, nbinsx, std::vector<double>(nbinsx+1, xmin).data(), nbinsy, std::vector<double>(nbinsy+1, ymin).data())
  {
    for (unsigned int ibin = 0; ibin <= nbinsx; ++ibin) {
      if (logx) {
        fBinXEdges.at(ibin) *= std::pow(xmax/xmin, double(ibin)/nbinsx);
      } else {
        fBinXEdges.at(ibin) += ibin * (xmax - xmin) / nbinsx;
      }
    } 

    for (unsigned int ibin = 0; ibin <= nbinsy; ++ibin) {
      if (logy) {
        fBinYEdges.at(ibin) *= std::pow(ymax/ymin, double(ibin)/nbinsy);
      } else {
        fBinYEdges.at(ibin) += ibin * (ymax - ymin) / nbinsy;
      }
    } 
  }
};

int
main(
  int argc,
  char** argv
)
{
  // disable printout messages from ROOT
  gErrorIgnoreLevel = kWarning;

  std::vector<std::string> fileNames(argv+1, argv+argc);
  auto data = ROOT::RDataFrame("data", fileNames)
    .Define("sixToDoubleChi2Dif", "ghSixParFit.chi2 - ghDoubleFit.chi2")
    .Define("deltaXmax", "std::fabs(ghDoubleFit.xmax_1-ghDoubleFit.xmax_2)")
    .Define("E", "std::pow(10, lgE)")
    .Define("dx0_1", "ghDoubleFit.xmax_1 - ghDoubleFit.x0_1")
    .Define("dx0_2", "ghDoubleFit.xmax_2 - ghDoubleFit.x0_2")
    .Define("dx0_a", "ghDoubleFit.w >= 0.5? dx0_1 : dx0_2")
    .Define("dx0_b", "ghDoubleFit.w >= 0.5? dx0_2 : dx0_1")
    .Define("lambda_a", "ghDoubleFit.w >= 0.5? ghDoubleFit.lambda_1 : ghDoubleFit.lambda_2")
    .Define("lambda_b", "ghDoubleFit.w >= 0.5? ghDoubleFit.lambda_2 : ghDoubleFit.lambda_1")
  ;

  // ? Data filters: select showers that are likely to have an anomalous profiles

  // - select showers whose both fitted maxima are within the fitted range
  auto cutXmaxInFov = data.Filter(
    "minDepth <= std::min(ghDoubleFit.xmax_1, ghDoubleFit.xmax_2) && " \
    "maxDepth >= std::max(ghDoubleFit.xmax_1, ghDoubleFit.xmax_2)"
  );

  // - selectshowers whose distance between maxima is above 250 g/cm^2
  auto cutDeltaXmax = cutXmaxInFov.Filter("deltaXmax > 250");
  

  // - select showers whose subshowers carry > 20% of Ecal each
  auto cutEcalMin = cutDeltaXmax.Filter(
    "0.2 <= ghDoubleFit.w && ghDoubleFit.w <= 0.8"
  );

  // - select showers based on shape parameters lambda_i and x0_i
  auto cutShape = cutEcalMin.Filter(
    "std::max(dx0_a, dx0_b) < 2e4 && " \
    "std::min(dx0_a, dx0_b) > 1e2 && " \
    "std::max(lambda_a, lambda_b) < 2e2 && " \
    "std::min(lambda_a, lambda_b) > 1"
  );

  // - select showers whose description though a double gaisser hillas
  // - is at least as good as the six-parameter gaisser-hillas
  auto cutChi2Improve = cutShape.Filter(
    "sixToDoubleChi2Dif > 1"
  );

  auto cutNinflec = data.Filter("ninflec > 2");

  // - wrap all cuts into a list
  using Filter_t = decltype(data.Filter(""));
  using FilterRef_t = std::reference_wrapper<Filter_t>;
  std::list<std::pair<std::string, FilterRef_t>> cutList;
  cutList.emplace_back("XmaxInFov", cutXmaxInFov);
  cutList.emplace_back("DeltaXmax", cutDeltaXmax);
  cutList.emplace_back("EcalMin", cutEcalMin);
  cutList.emplace_back("Shape", cutShape);
  cutList.emplace_back("Chi2Improve", cutChi2Improve);
  cutList.emplace_back("Ninflec", cutNinflec);

  // - reference to the final selection
  auto& cutFinal = cutChi2Improve;

  // ? Histograms: book the filling of histograms

  // - chi^2 distributions for the different fit functions
  TH1DLogModel th1Model_chi2("", "", 100, 1.5e-6, 3e3);

  auto histChi2_Double_NoCuts = data.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_XmaxInFov = cutXmaxInFov.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_DeltaXmax = cutDeltaXmax.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_EcalMin = cutEcalMin.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_Shape = cutShape.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_Chi2Improve = cutChi2Improve.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");
  auto histChi2_Double_Ninflec = cutNinflec.Histo1D(th1Model_chi2, "ghDoubleFit.chi2");

  auto histChi2_SixPar_NoCuts = data.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_XmaxInFov = cutXmaxInFov.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_DeltaXmax = cutDeltaXmax.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_EcalMin = cutEcalMin.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_Shape = cutShape.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_Chi2Improve = cutChi2Improve.Histo1D(th1Model_chi2, "ghSixParFit.chi2");
  auto histChi2_SixPar_Ninflec = cutNinflec.Histo1D(th1Model_chi2, "ghSixParFit.chi2");

  // - chi^2 difference between the six-parameter fit and the double-gh fit
  TH1DLogModel th1Model_chi2Dif("", "", 100, 3e-9, 2e3);

  auto histChi2_Dif_NoCuts = data.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_XmaxInFov = cutXmaxInFov.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_DeltaXmax = cutDeltaXmax.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_EcalMin = cutEcalMin.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_Shape = cutShape.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_Chi2Improve = cutChi2Improve.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");
  auto histChi2_Dif_Ninflec = cutNinflec.Histo1D(th1Model_chi2Dif, "sixToDoubleChi2Dif");

  // - x0 distributions after the Ecal and Xmax cuts
  TH1DLogModel th1Model_x0("", "", 100, 2e2, 1.2e6);
  auto histX0_a = cutEcalMin.Histo1D(th1Model_x0, "dx0_a");
  auto histX0_b = cutEcalMin.Histo1D(th1Model_x0, "dx0_b");

  // - lambda distributions after the Ecal and Xmax cuts
  TH1DLogModel th1Model_lambda("", "", 100, 4e-2, 2e4);
  auto histLambda_a = cutEcalMin.Histo1D(th1Model_lambda, "lambda_a");
  auto histLambda_b = cutEcalMin.Histo1D(th1Model_lambda, "lambda_b");

  // - distribution of lgE
  TH1DLogModel th1Model_lgE("", "", 20, 1e17, 1e20);
  auto histLgE_NoCuts = data.Histo1D(th1Model_lgE, "E");
  auto histLgE_Chi2Improve = cutChi2Improve.Histo1D(th1Model_lgE, "E");
  auto histLgE_Ninflec = cutNinflec.Histo1D(th1Model_lgE, "E");

  // ? Cut statistics (book and compute)

  std::list<std::pair<std::string, decltype(data.Count())>> filterCountList;
  filterCountList.emplace_back("NoCuts", data.Count());
  for (const auto& [name, cut] : cutList) {
    filterCountList.emplace_back(name, cut.get().Count());
  }

  for (auto& [name, count] : filterCountList) {
    const double n = count.GetValue();
    std::cout
      << std::setw(15) << name
      << std::setw(15) << n
      << std::setw(15) << n/filterCountList.front().second.GetValue()*100 << "%"
      << std::endl;
  }

  // ? Plot histograms

  MultiPageCanvas canvas("plots.pdf");

  // - chi2 histograms of the double Gaisser-Hillas fits

  // * Cut based on no cuts
  histChi2_Double_NoCuts->SetFillColorAlpha(kAzure+3, 0.3);
  histChi2_Double_NoCuts->SetLineColorAlpha(kBlack, 0);
  histChi2_Double_NoCuts->Draw("hist");
  // // * Cut based on Xmax in FoV
  // histChi2_Double_XmaxInFov->SetLineColorAlpha(kAzure-5, 0.8);
  // histChi2_Double_XmaxInFov->SetLineWidth(2);
  // histChi2_Double_XmaxInFov->Draw("hist same");
  // // * Cut based on distance between maxima
  // histChi2_Double_DeltaXmax->SetLineColorAlpha(kAzure-3, 0.8);
  // histChi2_Double_DeltaXmax->SetLineWidth(2);
  // histChi2_Double_DeltaXmax->Draw("hist same");
  // // * Cut based on minimum energy in each subshower
  // histChi2_Double_EcalMin->SetLineColorAlpha(kAzure-1, 0.8);
  // histChi2_Double_EcalMin->SetLineWidth(2);
  // histChi2_Double_EcalMin->Draw("hist same");
  // * Cut based on shape parameters
  histChi2_Double_Shape->SetLineColorAlpha(kAzure+3, 1.0);
  histChi2_Double_Shape->SetLineWidth(2);
  histChi2_Double_Shape->Draw("hist same");
  // * Cut based on minimum chi2 improvement
  histChi2_Double_Chi2Improve->SetLineColorAlpha(kAzure+3, 0.6);
  histChi2_Double_Chi2Improve->SetLineWidth(2);
  histChi2_Double_Chi2Improve->Draw("hist same");
  // * Cut based on the number of inflection points
  histChi2_Double_Ninflec->SetLineColorAlpha(kBlack, 0.7);
  histChi2_Double_Ninflec->SetLineWidth(2);
  histChi2_Double_Ninflec->Draw("hist same");

  canvas.Print(true, true);

  // - chi2 histograms of the six-parameter Gaisser-Hillas fits

  // * Double gaisser hillas: no cuts
  histChi2_SixPar_NoCuts->SetFillColorAlpha(kPink-7, 0.3);
  histChi2_SixPar_NoCuts->SetLineColorAlpha(kBlack, 0);
  histChi2_SixPar_NoCuts->Draw("hist");
  // // * Double gaisser hillas: Xmax in FoV
  // histChi2_SixPar_XmaxInFov->SetLineColorAlpha(kPink+2, 0.8);
  // histChi2_SixPar_XmaxInFov->SetLineWidth(2);
  // histChi2_SixPar_XmaxInFov->Draw("hist same");
  // // * Double gaisser hillas: distance between maxima
  // histChi2_SixPar_DeltaXmax->SetLineColorAlpha(kPink+3, 0.8);
  // histChi2_SixPar_DeltaXmax->SetLineWidth(2);
  // histChi2_SixPar_DeltaXmax->Draw("hist same");
  // // * Double gaisser hillas: minimum energy in each subshower
  // histChi2_SixPar_EcalMin->SetLineColorAlpha(kPink+4, 0.8);
  // histChi2_SixPar_EcalMin->SetLineWidth(2);
  // histChi2_SixPar_EcalMin->Draw("hist same");
  // * Cut based on shape parameters
  histChi2_SixPar_Shape->SetLineColorAlpha(kPink-7, 1.0);
  histChi2_SixPar_Shape->SetLineWidth(2);
  histChi2_SixPar_Shape->Draw("hist same");
  // * Cut based on minimum chi2 improvement
  histChi2_SixPar_Chi2Improve->SetLineColorAlpha(kPink-7, 0.8);
  histChi2_SixPar_Chi2Improve->SetLineWidth(2);
  histChi2_SixPar_Chi2Improve->Draw("hist same");
  // * Cut based on the number of inflection points
  histChi2_SixPar_Ninflec->SetLineColorAlpha(kBlack, 0.7);
  histChi2_SixPar_Ninflec->SetLineWidth(2);
  histChi2_SixPar_Ninflec->Draw("hist same");

  canvas.Print(true, true);

  // - chi2 histograms of both fits after the final cuts

  // * Double gaisser hillas
  histChi2_Double_NoCuts->Draw("hist");
  histChi2_Double_Shape->Draw("hist same");

  // * Six-parameter gaisser hillas
  histChi2_SixPar_NoCuts->SetFillStyle(kNone);
  histChi2_SixPar_NoCuts->SetLineColorAlpha(kPink-7, 0.7);
  histChi2_SixPar_NoCuts->SetLineWidth(3);
  histChi2_SixPar_NoCuts->Draw("hist same");
  histChi2_SixPar_Shape->Draw("hist same");

  canvas.Print(true, true);

  // - chi2 difference
  histChi2_Dif_NoCuts->SetLineColorAlpha(kBlack, 0.0);
  histChi2_Dif_NoCuts->SetFillColorAlpha(kOrange+8, 0.5);
  histChi2_Dif_NoCuts->Draw("hist");

  histChi2_Dif_Shape->SetLineColorAlpha(kOrange+10, 0.8);
  histChi2_Dif_Shape->SetLineWidth(2);
  histChi2_Dif_Shape->Draw("hist same");

  histChi2_Dif_Ninflec->SetLineColorAlpha(kBlack, 0.7);
  histChi2_Dif_Ninflec->SetLineWidth(2);
  histChi2_Dif_Ninflec->Draw("hist same");

  canvas.Print(true, true);

  // - parameters

  histX0_a->SetFillColorAlpha(kTeal+4, 0.6);
  histX0_a->SetLineColorAlpha(kBlack, 0.0);
  histX0_a->Draw("hist");

  histX0_b->SetLineColor(kViolet+2);
  histX0_b->SetLineWidth(2);
  histX0_b->Draw("hist same");

  canvas.Print(true, true);

  histLambda_a->SetFillColorAlpha(kTeal+4, 0.6);
  histLambda_a->SetLineColorAlpha(kBlack, 0.0);
  histLambda_a->Draw("hist");

  histLambda_b->SetLineColor(kViolet+2);
  histLambda_b->SetLineWidth(2);
  histLambda_b->Draw("hist same");

  canvas.Print(true, true);



  // - fraction of double anomalous showers versus energy
  histLgE_Chi2Improve->Divide(histLgE_NoCuts.GetPtr());
  histLgE_Chi2Improve->Scale(100);
  histLgE_Chi2Improve->SetLineWidth(2);
  histLgE_Chi2Improve->SetLineColor(kRed);
  histLgE_Chi2Improve->GetYaxis()->SetRangeUser(0.0, 1.05*histLgE_Chi2Improve->GetBinContent(1));
  histLgE_Chi2Improve->Draw("hist l");

  histLgE_Ninflec->Divide(histLgE_NoCuts.GetPtr());
  histLgE_Ninflec->Scale(100);
  histLgE_Ninflec->SetLineWidth(2);
  histLgE_Ninflec->SetLineColor(kBlue);
  histLgE_Ninflec->Draw("hist l same");

  canvas.Print(true, false);

  

  // return 0;

  // ? Plot anomalous profiles

  // auto profilePlotter = [&](
  //   // the fitted profile
  //   TGraph& dedx,
  //   // limits of the fitted range
  //   const double minDepth,
  //   const double maxDepth,
  //   // parameters of the double-gaisser-hillas fit
  //   const double dbl_ecal,
  //   const double dbl_w,
  //   const double dbl_x0_1,
  //   const double dbl_xmax_1,
  //   const double dbl_lambda_1,
  //   const double dbl_x0_2,
  //   const double dbl_xmax_2,
  //   const double dbl_lambda_2,
  //   // parameters of the six-param gaisser-hillas fit
  //   const double six_ymax,
  //   const double six_x0,
  //   const double six_xmax,
  //   const double six_p1,
  //   const double six_p2,
  //   const double six_p3
  // ) {
  //   dedx.SetPoint(dedx.GetN(), dedx.GetX()[dedx.GetN()-1], 0);
  //   dedx.SetPoint(dedx.GetN(), 0, 0);
  //   dedx.SetFillColorAlpha(kMagenta+4, 0.3);
  //   dedx.SetLineColorAlpha(kBlack, 0.0);
  //   dedx.Draw("af");

  //   GHDoubleFcn ghDouble;
  //   ghDouble.SetParameters(dbl_ecal, dbl_w, dbl_x0_1, dbl_xmax_1, dbl_lambda_1, dbl_x0_2, dbl_xmax_2, dbl_lambda_2);
  //   ghDouble.SetLineColorAlpha(kCyan+4, 0.7);
  //   ghDouble.SetLineWidth(3);
  //   ghDouble.SetLineStyle(kSolid);
  //   ghDouble.SetRange(minDepth, maxDepth);
  //   ghDouble.Draw("same l");

  //   GHSingleFcn ghFirstBump;
  //   ghFirstBump.SetParameters(dbl_ecal*dbl_w, dbl_x0_1, dbl_xmax_1, dbl_lambda_1);
  //   ghFirstBump.SetLineColorAlpha(kCyan+4, 0.7);
  //   ghFirstBump.SetLineWidth(2);
  //   ghFirstBump.SetLineStyle(kDashed);
  //   ghFirstBump.SetRange(minDepth, maxDepth);
  //   ghFirstBump.Draw("same l");

  //   GHSingleFcn ghSecondBump;
  //   ghSecondBump.SetParameters(dbl_ecal*(1-dbl_w), dbl_x0_2, dbl_xmax_2, dbl_lambda_2);
  //   ghSecondBump.SetLineColorAlpha(kCyan+4, 0.7);
  //   ghSecondBump.SetLineWidth(2);
  //   ghSecondBump.SetLineStyle(kDashed);
  //   ghSecondBump.SetRange(minDepth, maxDepth);
  //   ghSecondBump.Draw("same l");

  //   GHSixParamFcn ghSixPar;
  //   ghSixPar.SetParameters(six_ymax, six_x0, six_xmax, six_p1, six_p2, six_p3);
  //   ghSixPar.SetLineColorAlpha(kOrange+10, 0.7);
  //   ghSixPar.SetLineWidth(2);
  //   ghSixPar.SetLineStyle(kSolid);
  //   ghSixPar.SetRange(minDepth, maxDepth);
  //   ghSixPar.Draw("same l");

  //   canvas.Print();
  // };

  // ROOT::RDF::ColumnNames_t profilePlotterColumns = {
  //   // the fitted profile
  //   "dedx",
  //   // limits of the fitted range
  //   "minDepth",
  //   "maxDepth",
  //   // parameters of the double-gaisser-hillas fit
  //   "ghDoubleFit.ecal",
  //   "ghDoubleFit.w",
  //   "ghDoubleFit.x0_1",
  //   "ghDoubleFit.xmax_1",
  //   "ghDoubleFit.lambda_1",
  //   "ghDoubleFit.x0_2",
  //   "ghDoubleFit.xmax_2",
  //   "ghDoubleFit.lambda_2",
  //   // parameters of the six-param gaisser-hillas fit
  //   "ghSixParFit.ymax",
  //   "ghSixParFit.x0",
  //   "ghSixParFit.xmax",
  //   "ghSixParFit.p1",
  //   "ghSixParFit.p2",
  //   "ghSixParFit.p3"
  // };

  // cutChi2Improve.Foreach(profilePlotter, profilePlotterColumns);

  return 0;
}
