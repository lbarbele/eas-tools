#include <string>
#include <iostream>

#include "conex-file.h"
#include "conex-shower.h"

namespace conex {

  file::file(
    const std::string& fname
  ) :
    TFile(fname.c_str(), "read")
  {
    // check if file is open
    if (!is_open()) {
      std::cerr << "unable to open conex file " << fname << std::endl;
      return;
    }

    // get the trees and check
    m_header_tree.reset(Get<TTree>("Header"));
    m_shower_tree.reset(Get<TTree>("Shower"));

    if (!m_header_tree || !m_shower_tree) {
      Close();
      m_header_tree.reset(nullptr);
      m_shower_tree.reset(nullptr);
      return;
    }

    // set addresses of the header tree
    m_header_tree->SetBranchAddress("Seed1",&m_header.Seed1);
    m_header_tree->SetBranchAddress("Particle",&m_header.Particle);
    m_header_tree->SetBranchAddress("Alpha",&m_header.Alpha);
    m_header_tree->SetBranchAddress("lgEmin",&m_header.lgEmin);
    m_header_tree->SetBranchAddress("lgEmax",&m_header.lgEmax);
    m_header_tree->SetBranchAddress("zMin",&m_header.zMin);
    m_header_tree->SetBranchAddress("zMax",&m_header.zMax);
    m_header_tree->SetBranchAddress("Version",&m_header.Version);
    m_header_tree->SetBranchAddress("OutputVersion",&m_header.OutputVersion);
    m_header_tree->SetBranchAddress("HEModel",&m_header.HEModel);
    m_header_tree->SetBranchAddress("LEModel",&m_header.LEModel);
    m_header_tree->SetBranchAddress("HiLowEgy",&m_header.HiLowEgy);
    m_header_tree->SetBranchAddress("hadCut",&m_header.hadCut);
    m_header_tree->SetBranchAddress("emCut",&m_header.emCut);
    m_header_tree->SetBranchAddress("hadThr",&m_header.hadThr);
    m_header_tree->SetBranchAddress("muThr",&m_header.muThr);
    m_header_tree->SetBranchAddress("emThr",&m_header.emThr);
    m_header_tree->SetBranchAddress("haCut",&m_header.haCut);
    m_header_tree->SetBranchAddress("muCut",&m_header.muCut);
    m_header_tree->SetBranchAddress("elCut",&m_header.elCut);
    m_header_tree->SetBranchAddress("gaCut",&m_header.gaCut);

    m_header_tree->SetBranchAddress("lambdaLgE",m_header.lambdaLgE);
    m_header_tree->SetBranchAddress("lambdaProton",m_header.lambdaProton);
    m_header_tree->SetBranchAddress("lambdaPion",m_header.lambdaPion);
    m_header_tree->SetBranchAddress("lambdaHelium",m_header.lambdaHelium);
    m_header_tree->SetBranchAddress("lambdaNitrogen",m_header.lambdaNitrogen);
    m_header_tree->SetBranchAddress("lambdaIron",m_header.lambdaIron);

    if (m_header_tree->GetNbranches() > 27) {
      m_header_tree->SetBranchAddress("resamplingMode",&m_header.resamplingMode);
      m_header_tree->SetBranchAddress("modThreshold",&m_header.modThreshold);
      m_header_tree->SetBranchAddress("f19_cx",&m_header.f19_cx);
      m_header_tree->SetBranchAddress("f19_meson",&m_header.f19_meson);
      m_header_tree->SetBranchAddress("f19",&m_header.f19);
    }

    // set addresses of the shower tree
    m_shower_tree->SetBranchAddress("lgE",&m_shower.lgE);
    m_shower_tree->SetBranchAddress("zenith",&m_shower.zenith);
    m_shower_tree->SetBranchAddress("azimuth",&m_shower.azimuth);
    m_shower_tree->SetBranchAddress("Seed2",&m_shower.Seed2);
    m_shower_tree->SetBranchAddress("Seed3",&m_shower.Seed3);
    m_shower_tree->SetBranchAddress("Xfirst",&m_shower.Xfirst);
    m_shower_tree->SetBranchAddress("Hfirst",&m_shower.Hfirst);
    m_shower_tree->SetBranchAddress("XfirstIn",&m_shower.XfirstIn);
    m_shower_tree->SetBranchAddress("altitude",&m_shower.altitude);
    m_shower_tree->SetBranchAddress("X0",&m_shower.X0);
    m_shower_tree->SetBranchAddress("Xmax",&m_shower.Xmax);
    m_shower_tree->SetBranchAddress("Nmax",&m_shower.Nmax);
    m_shower_tree->SetBranchAddress("p1",&m_shower.p1);
    m_shower_tree->SetBranchAddress("p2",&m_shower.p2);
    m_shower_tree->SetBranchAddress("p3",&m_shower.p3);
    m_shower_tree->SetBranchAddress("chi2",&m_shower.chi2);
    m_shower_tree->SetBranchAddress("Xmx",&m_shower.Xmx);
    m_shower_tree->SetBranchAddress("Nmx",&m_shower.Nmx);
    m_shower_tree->SetBranchAddress("XmxdEdX",&m_shower.XmxdEdX);
    m_shower_tree->SetBranchAddress("dEdXmx",&m_shower.dEdXmx);
    m_shower_tree->SetBranchAddress("cpuTime",&m_shower.cpuTime);
    m_shower_tree->SetBranchAddress("nX",&m_shower.nX);

    m_shower_tree->SetBranchAddress("X",m_shower.X);
    m_shower_tree->SetBranchAddress("N",m_shower.N);
    m_shower_tree->SetBranchAddress("H",m_shower.H);
    m_shower_tree->SetBranchAddress("D",m_shower.D);
    m_shower_tree->SetBranchAddress("dEdX",m_shower.dEdX);
    m_shower_tree->SetBranchAddress("Mu",m_shower.Mu);
    m_shower_tree->SetBranchAddress("Gamma",m_shower.Gamma);
    m_shower_tree->SetBranchAddress("Electrons",m_shower.Electrons);
    m_shower_tree->SetBranchAddress("Hadrons",m_shower.Hadrons);
    m_shower_tree->SetBranchAddress("dMu",m_shower.dMu);
    m_shower_tree->SetBranchAddress("EGround",m_shower.EGround);

    // actually read the header
    m_header_tree->GetEntry(0);
  }

  bool
  file::is_open()
  const
  {
    return IsOpen() && !IsZombie() && !TestBit(kRecovered);
  }

  const header&
  file::get_header()
  const
  {
    return m_header;
  }

  long long
  file::get_n_showers()
  const
  {
    return m_shower_tree->GetEntries();
  }

  const shower&
  file::get_shower(
    size_t pos
  )
  {
    m_shower_tree->GetEntry(pos);
    return m_shower;
  }

  const shower&
  file::operator[](
    size_t pos
  )
  {
    return get_shower(pos);
  }

} // namespace conex