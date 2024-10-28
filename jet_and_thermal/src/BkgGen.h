#ifndef BkgGen_h
#define BkgGen_h

// Code slighted refactored from J. Putschke in msg, 2023.06.27

#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"
#include "TF1.h"
#include "TTree.h"


#include <string>
#include <vector>

// 2024.07.05  -- update to just read the background particles from 
// the TTree for the hydro simulation

struct BkgGen {
  double maxEta          { 1.1     };
  double dNdEta          { 650     };
  double T               { 0.291   };
  double minPtCut        { 0.2     };
  double PionMass        { 0.13957 };
  bool   include_neutral { true    }; // increases nParticles by ratio of 1.5
  unsigned int seed      { 0       };
  double chargedRatio    { 2./3.   };

  bool read_background_input = false;

  std::vector<int>   *bulk_PID {0};
  std::vector<float> *bulk_pt {0};
  std::vector<float> *bulk_eta {0};
  std::vector<float> *bulk_phi {0};
  std::vector<float> *bulk_E {0};

  TTree* tree;
  bool separate_bulk { false };
  bool is_init { false };
  int i_bulk_event = 0;
  int num_bulk_events = 0;
  void init(TTree* tree, bool _separate_bulk);

  // values to initialize from choices above
  // int nParticles { 0 }; // input: max_eta, minPtCut, T
  // TRandom3 rng   { };   // input: seed
  // TF1* fpt       { nullptr }; // <- new TF1("fpt","x*exp(-x/[0])",minPtCut,10.0);

  BkgGen() {};
  std::vector<fastjet::PseudoJet> operator()();
};

#endif
