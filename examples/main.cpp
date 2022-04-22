#include <iostream>
#include <iomanip>

#include <conex/file.h>
#include <corsika/longfile.h>

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2.h>

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

int main(void)
{
  const double minAge = 0.5;
  const double maxAge = 1.5;

  corsika::longfile file("/home/luan/Projetos/corsika-77410-extended-stackin/run/alphaEff/DAT000001.long");

  file.read();

  auto depths = file.get_particle_profile(corsika::profile_type::depth);
  auto charged = file.get_particle_profile(corsika::profile_type::charged);
  auto emioniz = file.get_deposit_profiles(corsika::profile_type::positron);
  auto muioniz = file.get_deposit_profiles(corsika::profile_type::muplus);
  auto haioniz = file.get_deposit_profiles(corsika::profile_type::hadron);

  TGraph g;

  const double xmax = file.get_fit().xmax;

  for (int i = 0; i < depths.size() - 1; ++i) {
    const double x = 0.5*(depths[i] + depths[i+1]);
    const double s = 3*x/(x + 2*xmax);
    if (s < minAge || s > maxAge) {
      continue;
    }
    const double dep = (emioniz[i] + muioniz[i] + haioniz[i]) / (depths[i+1] - depths[i]);
    const double n = charged[i];
    const double alpha = dep/n;
    g.SetPoint(g.GetN(), s, alpha);
  }

  TF1 fcnNerling("", alphaNerling, minAge, maxAge);
  TF1 fcnSong("", alphaSong, minAge, maxAge);

  TH2D hAxis("", "", 100, minAge, maxAge, 100, 0.0021, 0.00269999);
  hAxis.SetStats(false);

  TCanvas c;
  hAxis.Draw("axis");
  g.Draw("same l");
  fcnNerling.Draw("same l");
  fcnSong.Draw("same l");
  c.Print("alpha.pdf");

  return 0;
}

int
mainaa()
{
  const double minAge = 0.5;
  const double maxAge = 1.5;

  conex::file cxFile("../data/conex10302588_sibyll23d_629771534_100.root");

  TCanvas c;
  c.Print("alpha.pdf[");

  // plot ratio between song and nerling
  auto ratio = [](const double *s, const double *p){
    return alphaSong(s,p)/alphaNerling(s,p);
  };

  TF1 fcnRatio("", ratio, minAge, maxAge);
  fcnRatio.Draw("l");
  c.Print("alpha.pdf");

  TH2D hAxis("", "", 100, minAge, maxAge, 100, 0.0021, 0.0032999);
  hAxis.SetStats(false);

  TF1 fcnNerling("", alphaNerling, minAge, maxAge);
  TF1 fcnSong("", alphaSong, minAge, maxAge);

  for (int ishower = 0; ishower < cxFile.get_n_showers(); ++ishower) {
    const auto& shower = cxFile.get_shower(ishower);

    if (shower.get_lge() < 1117.1) {

      const double xmax = shower.get_xmax();

      const auto depths = shower.get_depths();
      const auto electrons = shower.get_electrons();
      const auto dedx = shower.get_dedx();

      TGraph g;

      for (int ipt = 0; ipt < shower.get_nx()-1; ++ipt) {
        const double x = 0.5*(depths[ipt] + depths[ipt+1]);
        const double s = 3*x/(x + 2*xmax);

        if (s < minAge || s > maxAge) {
          continue;
        }

        const double n = 0.5*(electrons[ipt] + electrons[ipt+1]);
        const double alpha = dedx[ipt] / n;
        
        g.SetPoint(g.GetN(), s, alpha);
      }

      hAxis.Draw("axis");
      g.Draw("same l");
      fcnNerling.Draw("same l");
      fcnSong.Draw("same l");
      c.Print("alpha.pdf");

      break;
    }
  }

  c.Print("alpha.pdf]");

  return 0;
}