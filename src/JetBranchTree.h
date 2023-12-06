#ifndef JetBranchTree__H
#define JetBranchTree__H

#include <TTree.h>
#include <TFile.h>
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>
#include "JTWalker.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

struct
JetBranch {
    //options
    bool fill_constituents { true };
    bool fill_area         { true };

    void reset();
    void fill(fastjet::PseudoJet& jet);

    //branches
    float pt;
    float phi;
    float eta;
    float area;
    int   charge;
    int   nconsts;

    // constituents
    vector<float> const_pt  {};
    vector<float> const_phi {};
    vector<float> const_eta {};
    vector<int>   const_charge {};
}



struct
JetTree {
  TTree* tree;

  // Control switches
  float jet_zg     {};
  float jet_Rg     {};
  float jet_mu     {};


  // Branches and branch options
  void init_branches();

  // the initiating quark and gluons
  float qg_pt  ; // 0 when not matched
  float qg_eta ; // 0 when not matched
  float qg_phi ; // 0 when not matched
  int   qg_id  ; // 0 when not matched

  //the truth jets  --  and missed
  bool   T_isQGmatched;
  bool   T_ismissed; //  not matched to Rjet
  float  Tjet_pt     ;
  float  Tjet_eta    ;
  float  Tjet_phi    ;
  int    Tjet_charge ;
  int    Tjet_nconst ;

  // truth get constituents
  vector<float> Tconst_pt  {};
  vector<float> Tconst_eta {};
  vector<float> Tconst_phi {};
  vector<int>   Tconst_charge {};

  //the reconstructed jets
  bool   R_isQGmatched ;  // is matched to a truth jet
  bool   R_isfake ; //  not matched to truth jet
  float  Rjet_pt     ;
  float  Rjet_eta    ;
  float  Rjet_phi    ;
  int    Rjet_charge ;
  int    Rjet_nconst ;

  // reconstructed get constituents
  vector<bool>  Rconst_isbkgd {}; // only knowable in simulation
  vector<float> Rconst_pt  {};
  vector<float> Rconst_eta {};
  vector<float> Rconst_phi {};
  vector<int>   Rconst_charge {};


  void add_matchedjet(fastjet::PseudoJet& qg_jet, fastjet::PseudoJet& reco_jet);

  void fill();

  void write();
  void clear_vectors();

  bool use_SD;
  const SD_criteria SD_crit; // in file JTWalker.h; only uses zcut, and beta,

  JetTree(string name="T", bool use_SD=false, SD_criteria sd_crit={});
};

#endif
