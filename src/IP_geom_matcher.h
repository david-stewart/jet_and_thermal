#ifndef IP_GEOM_MATCHER__H
#define IP_GEOM_MATCHER__H
#include "TH1D.h"
#include "TH2D.h"

// Code slighted refactored from J. Putschke in msg, 2023.06.27

// Read in lists of hadrons from Jetscape, match the leading IP
// and geometrically match to the leading jet

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include <vector>
#include "TTreeReader.h"
#include "JetClusterer.h"

class IP_geom_matcher {
    public:
    IP_geom_matcher(const std::string& fname);
    /* std::vector<fastjet::PseudoJet> operator()(); */

    bool has_matched_IP { false }; // looking at leading parton only
    fastjet::PseudoJet IP_jet {};
    fastjet::PseudoJet truth_jet {}; // matched jet to IP_jet
    std::vector<fastjet::PseudoJet> truth_jet_constituents {};
    std::vector<fastjet::PseudoJet> all_constituents {};
    fastjet::Selector select_circle = fastjet::SelectorCircle(0.4);

    float dR = -1.;

    const float IP_eta_max = 1.;
    const float IP_eta_min = -1.;

    bool print_QA = false;
    // void write();

    // std::array<double,3> deltas (fastjet::PseudoJet& a, fastjet::PseudoJet& b);
    // TH1D* LJ_deltaR { new TH1D("LJ_deltaR", "hg_deltaR", 60, 0., 6.) };
    // TH2D* LJ_delta_phieta { new TH2D("LJ_delta_phieta", ";#Delta#phi;#Delta#eta", 30, -M_PI, M_PI,  80, -4, 4.) };
    // TH1D* Edens_deltaR { new TH1D("Edens_deltaR", "Edens_deltaR;R;dpT/dphi/deta", 60., 0., 6.) };
    // TH2D* Edens_delta_phieta { new TH2D("Edens_delta_phieta", "pT density;#Delta#phi;#Delta#eta", 30, -M_PI, M_PI, 80, -4, 4.) };

// convenience functions
    fastjet::PseudoJet make_pseudojet(float pt, float eta, float phi, float E);

    bool return_false();

    TFile *fin;
    TTree* tree;

    long long int nEvents;
    long long int i_event {-1};

    bool next_event();
    bool cluster_match(JetClusterer& clusterer);

    // Declaration of leaf types
    Int_t           event;
    Double_t        sigma;
    Double_t        sigma_error;
    std::vector<int>     *PID;
    std::vector<float>   *pt;
    std::vector<float>   *eta;
    std::vector<float>   *phi;
    std::vector<float>   *E;
    std::vector<int>     *ip_PID;
    std::vector<float>   *ip_pt;
    std::vector<float>   *ip_eta;
    std::vector<float>   *ip_phi;
    std::vector<float>   *ip_E;

    // List of branches
    TBranch        *b_event;   //!
    TBranch        *b_sigma;   //!
    TBranch        *b_sigma_error;   //!
    TBranch        *b_PID;   //!
    TBranch        *b_pt;   //!
    TBranch        *b_eta;   //!
    TBranch        *b_phi;   //!
    TBranch        *b_E;   //!
    TBranch        *b_ip_PID;   //!
    TBranch        *b_ip_pt;   //!
    TBranch        *b_ip_eta;   //!
    TBranch        *b_ip_phi;   //!
    TBranch        *b_ip_E;   //!

    std::vector<fastjet::PseudoJet> operator()();
}; 


#endif
