#include <string>
#include <iostream>

#include <TFile.h>

#include <conex/file.h>
#include <conex/shower.h>

namespace conex {

  file::file(
    const std::vector<std::string>& fnames
  ) :
    m_shower_tree(std::make_unique<TChain>("Shower")),
    m_file_names(fnames)
  {
    // do nothing if the list of file names is empty
    if (fnames.empty()) {
      return;
    }

    // try to add the files to the TChain
    for (const auto& file_path : fnames) {
      if (m_shower_tree->Add(file_path.c_str(), 0) != 1) {
        std::cerr << "failed to open file " + file_path << std::endl;
        m_shower_tree.reset(nullptr);
        m_file_names.clear();
        return;
      }
    }

    // set the branch addresses of the Shower tree
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
  }
  
  bool
  file::is_open()
  const
  {
    return bool(m_shower_tree);
  }

  header
  file::get_header()
  const
  {
    // get pointer to the current TFile
    auto current_file = m_shower_tree->GetFile();

    // get the header tree from the current file
    auto header_tree = std::unique_ptr<TTree>(current_file->Get<TTree>("Header"));

    // create a header object
    header header_obj;

    // connect the header object with the header tree
    header_tree->SetBranchAddress("Seed1",&header_obj.Seed1);
    header_tree->SetBranchAddress("Particle",&header_obj.Particle);
    header_tree->SetBranchAddress("Alpha",&header_obj.Alpha);
    header_tree->SetBranchAddress("lgEmin",&header_obj.lgEmin);
    header_tree->SetBranchAddress("lgEmax",&header_obj.lgEmax);
    header_tree->SetBranchAddress("zMin",&header_obj.zMin);
    header_tree->SetBranchAddress("zMax",&header_obj.zMax);
    header_tree->SetBranchAddress("Version",&header_obj.Version);
    header_tree->SetBranchAddress("OutputVersion",&header_obj.OutputVersion);
    header_tree->SetBranchAddress("HEModel",&header_obj.HEModel);
    header_tree->SetBranchAddress("LEModel",&header_obj.LEModel);
    header_tree->SetBranchAddress("HiLowEgy",&header_obj.HiLowEgy);
    header_tree->SetBranchAddress("hadCut",&header_obj.hadCut);
    header_tree->SetBranchAddress("emCut",&header_obj.emCut);
    header_tree->SetBranchAddress("hadThr",&header_obj.hadThr);
    header_tree->SetBranchAddress("muThr",&header_obj.muThr);
    header_tree->SetBranchAddress("emThr",&header_obj.emThr);
    header_tree->SetBranchAddress("haCut",&header_obj.haCut);
    header_tree->SetBranchAddress("muCut",&header_obj.muCut);
    header_tree->SetBranchAddress("elCut",&header_obj.elCut);
    header_tree->SetBranchAddress("gaCut",&header_obj.gaCut);

    header_tree->SetBranchAddress("lambdaLgE",header_obj.lambdaLgE.data());
    header_tree->SetBranchAddress("lambdaProton",header_obj.lambdaProton.data());
    header_tree->SetBranchAddress("lambdaPion",header_obj.lambdaPion.data());
    header_tree->SetBranchAddress("lambdaHelium",header_obj.lambdaHelium.data());
    header_tree->SetBranchAddress("lambdaNitrogen",header_obj.lambdaNitrogen.data());
    header_tree->SetBranchAddress("lambdaIron",header_obj.lambdaIron.data());

    if (header_tree->GetNbranches() > 27) {
      header_tree->SetBranchAddress("resamplingMode",&header_obj.resamplingMode);
      header_tree->SetBranchAddress("modThreshold",&header_obj.modThreshold);
      header_tree->SetBranchAddress("f19_cx",&header_obj.f19_cx);
      header_tree->SetBranchAddress("f19_meson",&header_obj.f19_meson);
      header_tree->SetBranchAddress("f19",&header_obj.f19);

      header_obj.m_has_extensions = true;
    } else {
      header_obj.m_has_extensions = false;
    }

    // actually read the header
    header_tree->GetEntry(0);

    return header_obj;
  }

  long long
  file::get_n_showers()
  const
  {
    return m_shower_tree->GetEntries();
  }

  unsigned int
  file::get_n_files()
  const
  {
    return m_file_names.size();
  }

  int
  file::get_current_file_number()
  const
  {
    return is_open()? m_shower_tree->GetTreeNumber(): -1;
  }

  std::string
  file::get_current_file_name()
  const
  {
    return is_open()? m_shower_tree->GetFile()->GetName() : "";
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

  file::iterator
  file::begin()
  {
    return iterator(*this);
  }

  file::iterator
  file::end()
  {
    return iterator();
  }

} // namespace conex