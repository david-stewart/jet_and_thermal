#include "JetBranch.h"
#include <fastjet/contrib/RecursiveSymmetryCutBase.hh>  // for RecursiveSymm...
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>  // for FunctionOfPse...
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/RecursiveSymmetryCutBase.hh>  // for RecursiveSymm...
#include <fastjet/contrib/SoftDrop.hh>

using std::string;

void JetBranch::reset() {
    pt = 0.;
    phi =0.;
    eta =0.;
    if (fill_area) {
        area = 0.;
        ptlessarea = 0.;
    }
    if (fill_SD) {
        zg = 0;
        Rg = 0;
        mu = 0;
    }
    charge = 0.;
    nconsts = 0.;
    if (fill_constituents) {
        const_pt.clear();
        const_phi.clear();
        const_eta.clear();
        const_charge.clear();
        const_is_bkgd.clear();
    }
}

void JetBranch::fill(fastjet::PseudoJet& jet, float rho) {
    reset();
    pt = jet.perp();
    phi = jet.phi();
    eta = jet.eta();
    if (!jet.has_constituents()) {
        std::cout << " Error in src/JetBranch.cc: jet without constituents!" << std::endl;
        exit(1);
    }
    nconsts = 0;
    charge = 0;
    if (fill_area) {
        area = jet.area();
        ptlessarea = jet.perp()-area*rho;
    }
    for (auto& C : jet.constituents()) {
        if (fill_area && C.is_pure_ghost()) continue;
        nconsts += 1;
        auto index = (int) C.user_index();
        auto Q = index % 2;
        charge += Q;
        if (fill_constituents) {
            const_pt.push_back(C.perp());
            const_eta.push_back(C.eta());
            const_phi.push_back(C.phi());
            const_charge.push_back(Q);
            if (tag_bkgd) const_is_bkgd.push_back(fabs(index)>10);
        }
    }
}

void JetBranch::add_to_ttree(TTree* tree, std::string prefix, std::string postfix) {
    string pt_name = prefix + "pt" + postfix;
    string phi_name = prefix + "phi" + postfix;
    string eta_name = prefix + "eta" + postfix;
    string charge_name = prefix + "charge" + postfix;
    string nconst_name = prefix + "nconsts" + postfix;

    tree->Branch(pt_name.c_str(), &pt);
    tree->Branch(phi_name.c_str(), &phi);
    tree->Branch(eta_name.c_str(), &eta);
    if (fill_area) {
        string area_name = prefix + "area" + postfix;
        tree->Branch(area_name.c_str(), &area);

        string ptlessarea_name = prefix+"ptlessarea"+postfix;
        tree->Branch(ptlessarea_name.c_str(), &ptlessarea);
    }
    tree->Branch(charge_name.c_str(), &charge);
    tree->Branch(nconst_name.c_str(), &nconsts);

    if (fill_constituents) {
        string const_pt_name = prefix + "const_pt" + postfix;
        string const_phi_name = prefix + "const_phi" + postfix;
        string const_eta_name = prefix + "const_eta" + postfix;
        string const_charge_name = prefix + "const_charge" + postfix;

        tree->Branch(const_pt_name.c_str(), &const_pt);
        tree->Branch(const_phi_name.c_str(), &const_phi);
        tree->Branch(const_eta_name.c_str(), &const_eta);
        tree->Branch(const_charge_name.c_str(), &const_charge);

        if (tag_bkgd) {
            string name_bkgd = prefix + "const_bkgd" + postfix;
            tree->Branch(name_bkgd.c_str(), &const_is_bkgd);
        }
    }

    if (fill_SD) {
        string zg_name = prefix + "zg" + postfix;
        string Rg_name = prefix + "Rg" + postfix;
        string mu_name = prefix + "mu" + postfix;

        tree->Branch(zg_name.c_str(), &zg);
        tree->Branch(Rg_name.c_str(), &Rg);
        tree->Branch(mu_name.c_str(), &mu);
    }
}

