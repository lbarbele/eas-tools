#include <TGraph.h>
#include <TCanvas.h>

#include <corsika/longfile.h>

int
main(
  int argc,
  char** argv
)
{
  if (argc <= 1) {
    return 0;
  }

  TCanvas canvas;
  canvas.Print("corsika_profiles.pdf[");

  for (int ifile = 1; ifile < argc; ++ifile) {
    corsika::longfile longFile(argv[ifile]);

    for (auto profile : longFile) {
      auto depths = profile.get(corsika::profile::type::depth_dep);
      auto dedx = profile.get(corsika::profile::type::dedx_sum);

      TGraph g(depths.size(), &depths[0], &dedx[0]);
      g.Draw("al");
      g.SetLineWidth(2);
      g.SetLineColor(kBlue+4);
      canvas.Print("corsika_profiles.pdf");
    }
  }

  canvas.Print("corsika_profiles.pdf]");

  return 0;
}