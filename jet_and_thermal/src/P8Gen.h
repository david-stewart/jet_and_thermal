#ifndef P8Gen_h
#define P8Gen_h

// Code slighted refactored from J. Putschke in msg, 2023.06.27

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include <vector>
#include "TTreeReader.h"
#include "JetClusterer.h"

struct P8Gen {
  double sNN             { 200     };
  double pTHatMin        { 3.      };
  double pTHatMax        { 4.      };
  double maxEta          { 1.1     };
  double minPtCut        { 0.2     };
  double PionMass        { 0.13957 };
  bool   usePionMass     { false   }; // will use actual particle mass
  bool   collect_neutral { true    };
  bool   collect_charged { true    };
  std::string name_type  { "pp" };
  int    seed { 0 };
  int    nMaxBadGenerations = 10;

  double weightSum;
  double sigmaGen;

  Pythia8::Pythia pythia;

  P8Gen() {};
  bool is_init { false };
  void init();
  std::vector<fastjet::PseudoJet> operator()();

  std::vector<fastjet::PseudoJet> e5and6 {{}};
};

class P8TupReader {
    public:
    P8TupReader(const std::string& fname);
    std::vector<fastjet::PseudoJet> operator()();
    std::tuple<bool, fastjet::PseudoJet, std::vector<fastjet::PseudoJet>> lead_jet_only(JetClusterer&);
    bool next();

    TFile *fin;
    TTreeReader event_reader;
    TTreeReader part_reader;
    TTreeReaderValue<Float_t> px, py, pz, E, Xsec, XsecSigma, nPart;
}; 

#endif
