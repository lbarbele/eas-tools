#include <iostream>

#include <TCanvas.h>

#include <conex/file.h>

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
  canvas.Print("conex_profiles.pdf[");

  for (int ifile = 1; ifile < argc; ++ifile) {
    conex::file cxFile(argv[ifile]);

    for (const auto& shower : cxFile) {
      auto profile = shower.graph_dedx();
      profile.SetLineWidth(2);
      profile.Draw("al");
      canvas.Print("conex_profiles.pdf");
    }
  }

  canvas.Print("conex_profiles.pdf]");
}