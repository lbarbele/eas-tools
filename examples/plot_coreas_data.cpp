#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>

int
main(
  int argc,
  char** argv
)
{
  std::ifstream coreasBinsFile(argv[1]);

  std::string basePath = argv[1];
  if (basePath.find_last_of('/') != std::string::npos) {
    basePath = basePath.substr(0, basePath.find_last_of('/'));
  } else {
    basePath = ".";
  }

  double tFirst = 1e20;
  double tLast = 0;
  double eMin = 0;
  double eMax = 0;
  double eAbsMax = 0;

  while (coreasBinsFile.good()) {
    std::string antennaFileName;
    double x, y, z, zero, distance;
    coreasBinsFile >> antennaFileName >> x >> y >> z >> zero >> distance;
    if (!coreasBinsFile.good()) {
      break;
    }

    std::ifstream antennaFile(basePath + "/" + antennaFileName);

    double currentFirst = -1;
    double currentLast = 0;

    while(antennaFile.good()) {
      double t, ex, ey, ez, eabs;
      antennaFile >> t >> ex >> ey >> ez;
      if (!antennaFile.good()) {
        break;
      }
      eabs = std::sqrt(ex*ex + ey*ey + ez*ez);


      if (eabs != 0) {
        currentLast = t;
        if (currentFirst < 0) {
          currentFirst = t;
        }
      }

      if (eabs > eAbsMax) {
        eAbsMax = eabs;
      }

      if (std::max({ex, ey, ez}) > eMax) {
        eMax = std::max({ex, ey, ez});
      }

      if (std::min({ex, ey, ez}) < eMin) {
        eMin = std::min({ex, ey, ez});
      }
    }

    if (currentFirst < tFirst) {
      tFirst = currentFirst;
    }

    if (currentLast > tLast) {
      tLast = currentLast;
    }
  }

  coreasBinsFile.clear();
  coreasBinsFile.seekg(0);

  tLast = tFirst + 5e-6;

  gStyle->SetTitleFont(63, "XY");
  gStyle->SetLabelFont(63, "XY");
  gStyle->SetTitleSize(35, "XY");
  gStyle->SetLabelSize(25, "XY");
  gStyle->SetTitleYOffset(2);
  gStyle->SetTitleXOffset(1);
  gStyle->SetTitleSize(0.1, "T");

  TCanvas canvas("", "", 1600, 1600);
  canvas.SetRightMargin(0.03);
  canvas.SetBottomMargin(0.2);

  canvas.Divide(1,4, 0.01, 0.);
  canvas.Print("coreas.pdf[");

  while (coreasBinsFile.good()) {
    std::string antennaFileName;
    double x, y, z, zero, distance;
    coreasBinsFile >> antennaFileName >> x >> y >> z >> zero >> distance;
    if (!coreasBinsFile.good()) {
      break;
    }

    std::ifstream antennaFile(basePath + "/" + antennaFileName);

    TGraph gx, gy, gz, gabs;

    while (antennaFile.good()) {
      double t, ex, ey, ez, eabs;
      antennaFile >> t >> ex >> ey >> ez;
      if (!antennaFile.good()) {
        break;
      }

      if (t < tFirst) {
        continue;
      } else if (t > tLast) {
        break;
      }

      eabs = std::sqrt(ex*ex + ey*ey + ez*ez);

      gx.SetPoint(gx.GetN(), t*1e6, ex*1e9);
      gy.SetPoint(gy.GetN(), t*1e6, ey*1e9);
      gz.SetPoint(gz.GetN(), t*1e6, ez*1e9);
      gabs.SetPoint(gabs.GetN(), t*1e6, eabs*1e9);
    }

    canvas.cd(1);
    gx.Draw("al");
    gx.SetLineColor(kAzure+2);
    gx.GetYaxis()->CenterTitle();
    gx.GetYaxis()->SetTitle("E_{x} [nstatV/cm]");
    gx.GetYaxis()->SetRangeUser(eMin*1e9, eMax*1e9);

    canvas.cd(2);
    gy.Draw("al");
    gy.SetLineColor(kPink-3);
    gy.GetYaxis()->CenterTitle();
    gy.GetYaxis()->SetTitle("E_{y} [nstatV/cm]");
    gy.GetYaxis()->SetRangeUser(eMin*1e9, eMax*1e9);

    canvas.cd(3);
    gz.Draw("al");
    gz.SetLineColor(kSpring+5);
    gz.GetYaxis()->CenterTitle();
    gz.GetYaxis()->SetTitle("E_{z} [nstatV/cm]");
    gz.GetYaxis()->SetRangeUser(eMin*1e9, eMax*1e9);

    canvas.cd(4);
    gabs.Draw("al");
    gabs.GetYaxis()->SetRangeUser(0, eAbsMax*1e9);
    gabs.GetYaxis()->CenterTitle();
    gabs.GetYaxis()->SetTitle("|E| [nstatV/cm]");
    gabs.GetXaxis()->CenterTitle();
    gabs.GetXaxis()->SetTitle("t [#mus]");

    gx.SetTitle(antennaFileName.c_str());
    canvas.Print("coreas.pdf");
  }

  canvas.Print("coreas.pdf]");

  return 0;
}