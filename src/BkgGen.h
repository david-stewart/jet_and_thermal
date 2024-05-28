#ifndef BkgGen_h
#define BkgGen_h

// Code slighted refactored from J. Putschke in msg, 2023.06.27

#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"
#include "TF1.h"
#include <vector>

struct BkgGen {
  double maxEta          { 1.1     };
  double dNdEta          { 650     };
  double T               { 0.291   };
  double minPtCut        { 0.2     };
  double PionMass        { 0.13957 };
  bool   include_neutral { true    }; // increases nParticles by ratio of 1.5
  unsigned int seed      { 0       };
  double chargedRatio    { 2./3.   };

  bool is_init { false };
  void init();

  // values to initialize from choices above
  int nParticles { 0 }; // input: max_eta, minPtCut, T
  TRandom3 rng   { };   // input: seed
  TF1* fpt       { nullptr }; // <- new TF1("fpt","x*exp(-x/[0])",minPtCut,10.0);

  BkgGen() {};
  std::vector<fastjet::PseudoJet> operator()();
};

#endif
