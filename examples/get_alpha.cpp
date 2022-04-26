#include <iostream>
#include <iomanip>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>

#include <corsika/longfile.h>

double 
alphaNerling(
  const double* sptr,
  const double* p
)
{
  static const double c[5] = {3.90883, 1.05301, 9.91717, 2.41715, 0.13180};
  const double& s = *sptr;
  return (c[0]/std::pow(c[1] + s, c[2]) + c[3] + c[4]*s)/1000;
}

double
alphaSong(
  const double *s,
  const double *p
)
{
  // value obtained by Song et al. Astropart. Phys 14 (2000) 7
  // simulations with CORSIKA for gamma/proton/iron showers at 10^17 eV
  // they provide the values averaged over the entire shower development
  // return 0.002186; // gamma
  return 0.002193; // proton
  // return 0.002189; // iron
}

int main(int argc, char** argv)
{
  const double minAge = 0.5;
  const double maxAge = 1.5;

  TCanvas canvas;
  // canvas.Print("alpha.pdf[");

  TH2D hAxis("", "", 100, 0.45, 1.55, 100, 0.0021, 0.0028);
  hAxis.SetStats(false);
  hAxis.Draw("axis");

  for (int ifile = 1; ifile < argc; ++ifile) {

    corsika::longfile longFile(argv[ifile]);

    for (auto profile : longFile) {
      const int nx = profile.size();

      // shower xmax from the gaisser-hillas fit
      const double xmax = profile.get_fit().get_xmax();

      // slant depth at the middle of each longitudinal step
      auto depths = profile.get(corsika::profile::type::depth_dep);

      // shower age values at each of the above depths
      std::valarray<double> ages = 3/(1 + 2*xmax/depths);

      // particle count at the middle of each longitudinal step
      auto part = std::valarray<double>(nx);
      // part += profile.get(corsika::profile::type::electrons);
      // part += profile.get(corsika::profile::type::positrons);
      part += profile.get(corsika::profile::type::charged);
      part = 0.5 * (part + part.shift(-1));

      // energy deposit at the middle of each longitudinal step
      auto dedx = std::valarray<double>(nx);
      // dedx += profile.get(corsika::profile::type::gamma_dep);
      dedx += profile.get(corsika::profile::type::em_ioniz);
      // dedx += profile.get(corsika::profile::type::em_cut);
      dedx += profile.get(corsika::profile::type::mu_ioniz);
      // dedx += profile.get(corsika::profile::type::mu_cut);
      dedx += profile.get(corsika::profile::type::hadr_ioniz);
      // dedx += profile.get(corsika::profile::type::hadr_cut);
      // dedx += profile.get(corsika::profile::type::netrino_dep);
      // dedx += profile.get(corsika::profile::type::dedx_sum);
      dedx /= profile.get_step_size();

      // energy deposit per particle
      std::valarray<double> alpha = dedx/part;

      // find first/last elements based on age limits
      int first = 0;
      int last = nx-1;
      for (int i = 0; i < nx-1; ++i) {
        if (ages[i] < minAge) {
          ++first;
        } else {
          break;
        }
      }

      for (int i = nx-2; i >= 0; --i) {
        if (ages[i] > maxAge) {
          --last;
        } else {
          break;
        }
      }

      auto graph = new TGraph(last - first, &ages[first], &alpha[first]);
      graph->SetLineColorAlpha(kBlack, 0.05);
      graph->Draw("same l");
    }

  }

  TF1 fcnNerling("", alphaNerling, minAge, maxAge);
  fcnNerling.SetLineColor(kRed);
  fcnNerling.Draw("same l");

  TF1 fcnSong("", alphaSong, minAge, maxAge);
  fcnSong.SetLineColor(kBlue);
  fcnSong.Draw("same l");

  canvas.Print("alpha.pdf");

  return 0;
}