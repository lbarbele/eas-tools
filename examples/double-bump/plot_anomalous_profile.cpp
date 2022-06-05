#include <iostream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLegend.h>

#include <conex/file.h>
#include <corsika/longfile.h>
#include <corsika/binaryfile.h>

int
main(
  int argc,
  char** argv
)
{
  if (argc < 3) {
    return 1;
  }

  conex::file cnxFile(argv[1]);
  corsika::longfile cskLongFile(argv[2]);

  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gStyle->SetTitleFont(62, "XY");

  TCanvas canvas("", "", 900, 900);
  canvas.Print("anomalous.pdf[");

  // Draw the conex profile and the corsika profile together
  {
    auto cnxProfile = cnxFile.get_shower(cnxFile.get_n_showers() - 1).graph_dedx();

    auto cskData = *cskLongFile.begin();
    auto cskDepths = cskData.get(corsika::profile::type::depth_dep);
    auto cskdEdX = cskData.get(corsika::profile::type::dedx_sum);
    cskdEdX /= cskData.get_step_size();
    TGraph cskProfile(cskData.size()-3, &cskDepths[0], &cskdEdX[0]);

    cnxProfile.SetTitle("Longitudinal energy deposit profiles");
    cnxProfile.GetYaxis()->CenterTitle();
    cnxProfile.GetXaxis()->CenterTitle();
    cnxProfile.GetXaxis()->SetTitle("Slant depth [g cm^{-2}]");
    cnxProfile.GetYaxis()->SetTitle("dE/dX [GeV cm^{2} g^{-1}]");
    cnxProfile.SetLineColor(kAzure+2);
    cnxProfile.SetFillColorAlpha(kAzure+2, 0.4);
    cnxProfile.SetLineWidth(3);
    cnxProfile.Draw("afl");

    cskProfile.SetLineColorAlpha(kPink-3,0.75);
    cskProfile.SetLineWidth(3);
    cskProfile.Draw("same l");

    TLegend legend(0.65, 0.73, 0.83, 0.83);
    legend.AddEntry(&cnxProfile, "CONEX", "lf");
    legend.AddEntry(&cskProfile, "CORSIKA", "l");
    legend.Draw();

    canvas.Print("anomalous.pdf");
  }

  // draw the lateral distribution of particles in corsika
  if (argc > 3) {
    corsika::binaryfile cskBinaryFile(argv[3]);

    auto shower = *cskBinaryFile.begin();

    // find the core position
    double xCore = 0;
    double yCore = 0;
    double wpSum = 0;
    double wSum = 0;
    double tMean = 0;
    double tSigma = 0;
    for (const auto& particle : shower) {
      const int id = int(particle[0])/1000;
      const float& px = particle[1];
      const float& py = particle[2];
      const float& pz = particle[3];
      const float x = particle[4]/100000.;
      const float y = particle[5]/100000.;
      const float& t = particle[6];
      const float& w = particle[7];

      if (id <= 0 || id > 5656) {
        continue;
      }

      const double p = std::sqrt(px*px + py*py + pz*pz);

      xCore += w*p*x;
      yCore += w*p*y;
      wpSum += w*p;

      tMean += w*t;
      tSigma += w*t*t;
      wSum += w;
    }

    xCore /= wpSum;
    yCore /= wpSum;
    tMean /= wSum;
    tSigma = std::sqrt(tSigma/wSum - tMean*tMean);

    const double tBinSize = 250;

    yCore = 0;
    xCore = 0;
    TH2D lateral("", "", 100, xCore-3.999, xCore+3.999, 100, yCore-3.999, yCore+3.999);
    TH1D timeStructure("", "", 6*tSigma/tBinSize, tMean-3*tSigma, tMean+3*tSigma);
    TH1D timeElectrons("", "", 6*tSigma/tBinSize, tMean-3*tSigma, tMean+3*tSigma);
    TH1D timeMuons("", "", 6*tSigma/tBinSize, tMean-3*tSigma, tMean+3*tSigma);
    TH1D timeOther("", "", 6*tSigma/tBinSize, tMean-3*tSigma, tMean+3*tSigma);

    for (const auto& particle : shower) {
      const int id = int(particle[0])/1000;
      const float x = particle[4]/100000.;
      const float y = particle[5]/100000.;
      const float& t = particle[6];
      const float& w = particle[7];

      if (id <= 0 || id > 5656) {
        continue;
      }

      lateral.Fill(x, y, w);

      if (id == 1) {
        continue;
      }

      timeStructure.Fill(t, w);

      if (id == 2 || id == 3) {
        timeElectrons.Fill(t, w);
      } else if (id == 5 || id == 6) {
        timeMuons.Fill(t,w);
      } else {
        timeOther.Fill(t,w);
      }
    }

    gStyle->SetPalette(kSunset);
    lateral.SetStats(false);
    lateral.SetTitle("Particle density at ground");
    lateral.GetXaxis()->CenterTitle();
    lateral.GetYaxis()->CenterTitle();
    lateral.GetXaxis()->SetTitle("x [km]");
    lateral.GetYaxis()->SetTitle("y [km]");
    lateral.Draw("colz");
    canvas.SetLogz(true);
    canvas.Print("anomalous.pdf");
    canvas.SetLogz(false);

    timeStructure.SetStats(false);
    timeStructure.SetTitle("Distribution of arrival times");
    timeStructure.GetXaxis()->CenterTitle();
    timeStructure.GetYaxis()->CenterTitle();
    timeStructure.GetXaxis()->SetTitle("t [ns]");
    timeStructure.GetYaxis()->SetTitle("Particle count");
    timeStructure.SetLineColor(kBlue+4);
    timeStructure.SetLineWidth(2);
    timeStructure.Draw("hist");

    timeElectrons.SetLineColor(kAzure+2);
    timeElectrons.SetLineWidth(2);
    timeElectrons.Draw("hist same");

    timeMuons.SetLineColor(kPink-3);
    timeMuons.SetLineWidth(2);
    timeMuons.Draw("hist same");

    // timeOther.SetLineColor(kSpring+5);
    // timeOther.SetLineWidth(2);
    // timeOther.Draw("hist same");

    TLegend legend(0.65, 0.73, 0.83, 0.83);
    legend.AddEntry(&timeStructure, "Charged", "l");
    legend.AddEntry(&timeElectrons, "Electrons", "l");
    legend.AddEntry(&timeMuons, "Muons", "l");
    legend.Draw();

    canvas.Print("anomalous.pdf");
  }

  canvas.Print("anomalous.pdf]");

  return 0;
}