#ifndef JetBranchTree__H
#define JetBranchTree__H

#include <TTree.h>
#include <TFile.h>
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>
/* #include "JTWalker.h" */

#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"

using std::string;
using std::vector;

// Note the constituent indices are 0, 1, -1; for background constituents they are 100, 101, -101
enum JETBRANCH_OPTIONS {
    CONSTITUENTS, AREA, SOFT_DROP, TAG_BKGD, ANGULARITY, VEC_CONSTITUENTS
};

struct JetConstBranch {
  float dR(fastjet::PseudoJet& c, fastjet::PseudoJet& j);
  bool is_truth;
  vector<bool> c_is_matched;
  vector<int> cindex;
  vector<float> cpt{}, ceta{}, cphi{}, cdR{};
  float jpt{0.}, jeta{0.}, jphi{0.};
  // float area {0.};
  // float ptlessarea {0.};

  JetConstBranch(string tag, bool _is_truth, TTree* tree);
  void add_jet(fastjet::PseudoJet& jet, float rho=0.);
  void fill_matched(fastjet::PseudoJet& truth_jet, fastjet::PseudoJet& reco_jet);
};

struct
JetBranch {
    //options
    bool fill_constituents { false };
    bool fill_area         { false };
    bool fill_SD { false };
    bool tag_bkgd { false };
    bool fill_angularity { false };
    bool fill_vec_constituents { false };
    int  nbranch_const = 0;
    vector<float> vec_lead_const_pt {};

    JetBranch(const vector<JETBRANCH_OPTIONS> _={}, float _jet_R=0.4 ); 
    const float jet_R;
    /* bool _fill_constituents=true, bool _fill_area=true, bool */
    /*         _fill_SD=false, bool _tag_bkgd=false) */ 
    /*     : fill_constituents {_fill_constituents} */
    /*     , fill_area{_fill_area} */
    /*     , fill_SD{_fill_SD} */
    /*     , tag_bkgd{_tag_bkgd} */ 
    /* {}; */

    void reset();
    void fill(fastjet::PseudoJet& jet, float& rho, bool is_firstjet, float xsec=1.0);
    void add_to_ttree(TTree* tree, std::string prefix, std::string postfix="");

    //branches
    float pt;
    float phi;
    float eta;
    float area;
    float angularity;
    // float ptlessarea;
    int   charge;
    int   nconsts;
    bool  isleadjet;

    // constituents
    vector<float> const_pt  {};
    vector<float> const_phi {};
    vector<float> const_eta {};
    vector<int>   const_charge {};
    vector<bool>  const_is_bkgd {};
    vector<float> vec_cpt{}; // for fragmentation

    //softdrop
    float zg;
    float Rg;
    float mu;

    // for splitting function
    TH1D* jet_w;
    TH2D* jet_const;
};

#endif
