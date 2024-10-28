#ifndef JetTree__H
#define JetTree__H

#include <TTree.h>
#include <TFile.h>
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>
/* #include "JTWalker.h" */

#include <string>
#include <vector>

using std::string;
using std::vector;

struct SD_criteria {
  float zg      ;
  float beta    ;
  float R0      ;
  float min_pt  ;
  SD_criteria ( float _zg=0.3, float _beta=0., float _R0=0.4, float _min_pt=1. ) :
    zg{_zg}
  , beta{_beta}
  , R0 {_R0}
  , min_pt { _min_pt}
  {}
};

struct
JetTree {
  TTree* tree;

  void init_branches();

  vector<float> qg_pt      {};
  vector<float> qg_eta     {};
  vector<float> qg_phi     {};
  vector<int>   qg_id      {};

  vector<float> jet_pt     {};
  vector<float> jet_eta    {};
  vector<float> jet_phi    {};
  vector<int>   jet_charge {};
  vector<int>   jet_nconst {};

  vector<float> jet_zg     {};
  vector<float> jet_Rg     {};
  vector<float> jet_mu     {};

  void add_jetpair(fastjet::PseudoJet& qg_jet, fastjet::PseudoJet& reco_jet);
  void fill();

  void write();
  void clear_vectors();

  bool use_SD;
  const SD_criteria SD_crit; // in file JTWalker.h; only uses zcut, and beta,

  JetTree(string name="T", bool use_SD=false, SD_criteria sd_crit={});
};

#endif
