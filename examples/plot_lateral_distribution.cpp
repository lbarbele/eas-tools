#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include <corsika/binaryfile.h>
#include <corsika/shower.h>

#include <util/constants.h>
#include <util/math.h>
#include <util/point.h>
#include <util/units.h>
#include <util/vector.h>

#include <TCanvas.h>
#include <TExec.h>
#include <TH2.h>
#include <TLatex.h>
#include <TStyle.h>

int main(
  int argc,
  char** argv
)
{
  using namespace units::literals;

  // * configs
  static const int nbins = 100;
  static const auto max_dist = 5_km;

  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gStyle->SetTitleFont(62, "XY");


  // * check command line
  if (argc <= 1) {
    std::cerr << "missing input files" << std::endl;
    return 0;
  }

  // * the output canvas

  TCanvas canvas("", "", 1000, 1000);
  canvas.Divide(2, 2);
  canvas.Print("lateral_profiles.pdf[");

  // * loop over input files
  for (int ifile = 1; ifile < argc; ++ifile) {
    // open the corsika file
    corsika::binaryfile csk_file(argv[ifile]);

    // shower counter (for printout)
    unsigned ishower = 0;


    // loop over showers in the current file
    for (const auto& shower : csk_file) {

      // increment the shower counter
      ++ishower;

      // create histograms for this shower
      TH2D gammas2d("", "gammas", nbins, -max_dist/1_km, max_dist/1_km, nbins, -max_dist/1_km, max_dist/1_km);
      TH2D electrons2d("", "electrons", 100, -max_dist/1_km, max_dist/1_km, 100, -max_dist/1_km, max_dist/1_km);
      TH2D muons2d("", "muons", nbins, -max_dist/1_km, max_dist/1_km, nbins, -max_dist/1_km, max_dist/1_km);
      TH2D others2d("", "other", nbins, -max_dist/1_km, max_dist/1_km, nbins, -max_dist/1_km, max_dist/1_km);

      // loop over particles in the current file
      for (const auto& particle : shower) {

        // - get particle info

        // particle id
        const auto id_data = std::lround(particle[0]);
        const auto id = id_data/1000;
        const auto hadronic_generation = (id_data%1000)/10;
        const auto i_observation_level = id_data%10;

        // particle momentum
        const util::vector_t<units::gev_per_c_t<double>> momentum(
          particle[0] * (1_GeV/1_c),
          particle[1] * (1_GeV/1_c),
          particle[2] * (1_GeV/1_c),
          util::frame::corsika_observer);

        // particle position > convert from curved to euclidean coordinates, if necessary
        util::point_t<units::centimeter_t<double>> position(
          particle[4] * 1_cm,
          particle[5] * 1_cm,
          shower.get_header()[46+i_observation_level] * 1_cm,
          util::frame::corsika_observer);

        // arrival time
        const units::nanosecond_t<double> time = particle[6] * 1_ns;

        // statistical weight
        const double weight = particle[7];

        // - process particle

        switch (id) {
          case 1: // gammas
            gammas2d.Fill(position[0]/1_km, position[1]/1_km, weight);
            break;
          case 2: // positrons
          case 3: // electrons
            electrons2d.Fill(position[0]/1_km, position[1]/1_km, weight);
            break;
          case 5: // muon+
          case 6: // muon-
            muons2d.Fill(position[0]/1_km, position[1]/1_km, weight);
            break;
          default:
            others2d.Fill(position[0]/1_km, position[1]/1_km, weight);
            break;
        }

      } // particles

      // * draw the histograms

      canvas.cd(1);
      gammas2d.SetStats(false);
      gammas2d.GetXaxis()->SetTitle("x [km]");
      gammas2d.GetYaxis()->SetTitle("y [km]");
      gammas2d.GetXaxis()->CenterTitle();
      gammas2d.GetYaxis()->CenterTitle();
      TExec exec1("exec1", "gStyle->SetPalette(kDeepSea);");
      gammas2d.Draw("axis");
      exec1.Draw();
      gammas2d.Draw("colz same");
      gPad->SetLogz(true);

      canvas.cd(2);
      electrons2d.SetStats(false);
      electrons2d.GetXaxis()->SetTitle("x [km]");
      electrons2d.GetYaxis()->SetTitle("y [km]");
      electrons2d.GetXaxis()->CenterTitle();
      electrons2d.GetYaxis()->CenterTitle();
      TExec exec2("exec2", "gStyle->SetPalette(kAvocado);");
      electrons2d.Draw("axis");
      exec2.Draw();
      electrons2d.Draw("colz same");
      gPad->SetLogz(true);

      canvas.cd(3);
      muons2d.SetStats(false);
      muons2d.GetXaxis()->SetTitle("x [km]");
      muons2d.GetYaxis()->SetTitle("y [km]");
      muons2d.GetXaxis()->CenterTitle();
      muons2d.GetYaxis()->CenterTitle();
      TExec exec3("exec3", "gStyle->SetPalette(kSunset);");
      muons2d.Draw("axis");
      exec3.Draw();
      muons2d.Draw("colz same");
      gPad->SetLogz(true);

      // TLatex description(0, 0, "teste");
      std::string description = argv[ifile] + std::string(":") + std::to_string(ishower);
      TLatex().DrawLatexNDC(0.002, 0.002, description.c_str());

      canvas.cd(4);
      others2d.SetStats(false);
      others2d.GetXaxis()->SetTitle("x [km]");
      others2d.GetYaxis()->SetTitle("y [km]");
      others2d.GetXaxis()->CenterTitle();
      others2d.GetYaxis()->CenterTitle();
      TExec exec4("exec4", "gStyle->SetPalette(kLightTemperature);");
      others2d.Draw("axis");
      exec4.Draw();
      others2d.Draw("colz same");
      gPad->SetLogz(true);

      canvas.Print("lateral_profiles.pdf");

    } // showers
  } // files

  canvas.Print("lateral_profiles.pdf]");
}