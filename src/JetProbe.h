#ifndef JetProbe__H
#define JetProbe__H

#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include "fastjet/PseudoJet.hh"
using fastjet::PseudoJet;

#include <string>
#include <vector>

using std::string;
using std::vector;
using std::pair;

struct
JetProbe {
    TRandom3& rnd3;
    pair<float,float> eta_range;

    //options
    int  nbranch_const = 0;
    JetProbe(TRandom3& _rand, float eta=0.7, float eta_2=666., float _pt=30.) 
        : rnd3 {_rand}
        , eta_range {-eta, eta}
        , pt {_pt}
    {
        if (eta_2 != 666.) eta_range.second = eta_2;
    };
    void reset();
    void fill(fastjet::PseudoJet& jet, float rho=0);
    void add_to_ttree(TTree* tree, std::string postfix="");

    vector<PseudoJet> operator()();
    PseudoJet jet {};

    //branches
    float pt;
    float phi;
    float eta;
    int   nconsts {1};

    // constituents
    vector<float> const_pt  {};
    vector<float> const_phi {};
    vector<float> const_eta {};
};

#endif
