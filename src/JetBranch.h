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

// Note the constituent indices are 0, 1, -1; for background constituents they are 100, 101, -101
enum OPTIONS {
    CONSTITUENTS, AREA, SOFT_DROP, TAG_BKGD, ANGULARITY
};

struct
JetBranch {
    //options
    bool fill_constituents { false };
    bool fill_area         { false };
    bool fill_SD { false };
    bool tag_bkgd { false };
    bool fill_angularity { false };

    JetBranch(const vector<OPTIONS> _={}, float _jet_R=0.4 ); 
    const float jet_R;
    /* bool _fill_constituents=true, bool _fill_area=true, bool */
    /*         _fill_SD=false, bool _tag_bkgd=false) */ 
    /*     : fill_constituents {_fill_constituents} */
    /*     , fill_area{_fill_area} */
    /*     , fill_SD{_fill_SD} */
    /*     , tag_bkgd{_tag_bkgd} */ 
    /* {}; */

    void reset();
    void fill(fastjet::PseudoJet& jet, float rho=0);
    void add_to_ttree(TTree* tree, std::string prefix, std::string postfix="");

    //branches
    float pt;
    float phi;
    float eta;
    float area;
    float angularity;
    float ptlessarea;
    int   charge;
    int   nconsts;

    // constituents
    vector<float> const_pt  {};
    vector<float> const_phi {};
    vector<float> const_eta {};
    vector<int>   const_charge {};
    vector<bool>  const_is_bkgd {};

    //softdrop
    float zg;
    float Rg;
    float mu;
};

#endif
