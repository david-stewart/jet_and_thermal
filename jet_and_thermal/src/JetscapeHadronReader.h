#ifndef JETSCAPEHADRONREADER__h
#define JETSCAPEHADRONREADER__h
#include "TH1D.h"
#include "TH2D.h"

// Code slighted refactored from J. Putschke in msg, 2023.06.27

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include <vector>
#include "TTreeReader.h"
#include "JetClusterer.h"

class JetscapeHadronReader {
    public:
    JetscapeHadronReader(const std::string& fname);
    /* std::vector<fastjet::PseudoJet> operator()(); */

    bool has_leading_jet { false };
    int  IP_match_index { -1 }; // -1 for no match, 0 for leading, 1 for subleading
    fastjet::PseudoJet IP_leading {};
    fastjet::PseudoJet IP_sub {};
    fastjet::PseudoJet leading_jet {};
    std::vector<fastjet::PseudoJet> leading_jet_constituents {};

    bool print_QA = false;
    void write();

    std::array<double,3> deltas (fastjet::PseudoJet& a, fastjet::PseudoJet& b);

    TH1D* LJ_deltaR { new TH1D("LJ_deltaR", "hg_deltaR", 60, 0., 6.) };
    TH2D* LJ_delta_phieta { new TH2D("LJ_delta_phieta", ";#Delta#phi;#Delta#eta", 30, -M_PI, M_PI,  80, -4, 4.) };
    TH1D* Edens_deltaR { new TH1D("Edens_deltaR", "Edens_deltaR;R;dpT/dphi/deta", 60., 0., 6.) };
    TH2D* Edens_delta_phieta { new TH2D("Edens_delta_phieta", "pT density;#Delta#phi;#Delta#eta", 30, -M_PI, M_PI, 80, -4, 4.) };

// convenience functions
    fastjet::PseudoJet make_pseudojet(float pt, float eta, float phi, float E);

    TFile *fin;
    TTree* tree;
    long long int nEvents;
    long long int i_event {-1};

    bool next_event();

    bool cluster_jets(JetClusterer& clusterer);

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
