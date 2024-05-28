#include "IP_geom_matcher.h"
#include <iostream>
#include "TRandom3.h"
#include "Pythia8/Info.h"
#include "TLorentzVector.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

// void IP_geom_matcher::write() {
    // LJ_deltaR ->Write();
    // LJ_delta_phieta ->Write();
    // Edens_deltaR ->Write();
    // Edens_delta_phieta ->Write();
// 
// }

IP_geom_matcher::IP_geom_matcher(const std::string& fname) 
    : fin { new TFile(fname.c_str(), "read") }
    , tree { (TTree*)fin->Get("T") }
    , nEvents { tree->GetEntriesFast() }
    , event {0}
    , sigma {0}
    , sigma_error {0}
    , PID {0}
    , pt {0}
    , eta {0}
    , phi {0}
    , E {0}
    , ip_PID {0}
    , ip_pt {0}
    , ip_eta {0}
    , ip_phi {0}
    , ip_E {0}
{
    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("sigma", &sigma);
    tree->SetBranchAddress("sigma_error", &sigma_error);
    tree->SetBranchAddress("PID", &PID);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("E", &E);
    tree->SetBranchAddress("ip_PID", &ip_PID);
    tree->SetBranchAddress("ip_pt", &ip_pt);
    tree->SetBranchAddress("ip_eta", &ip_eta);
    tree->SetBranchAddress("ip_phi", &ip_phi);
    tree->SetBranchAddress("ip_E", &ip_E);
} 

bool IP_geom_matcher::next_event() {
    i_event++;
    if (i_event >= nEvents) {
        return false;
    } 
    if (print_QA) cout << " FIXME Event: " << i_event << endl;
    return tree->GetEntry(i_event);
}

fastjet::PseudoJet IP_geom_matcher::make_pseudojet(float pt, float eta, float phi, float E) {
    TLorentzVector C;
    C.SetPtEtaPhiE(pt, eta, phi, E);
    return fastjet::PseudoJet(C);
}

bool IP_geom_matcher::return_false() {
    has_matched_IP = false;
    truth_jet = fastjet::PseudoJet();
    truth_jet_constituents.clear();
    return false;
}

bool IP_geom_matcher::cluster_match(JetClusterer& clusterer) {
    int n_IP = ip_pt->size();
    if (n_IP < 1) {
        if (print_QA) cout << "FIXME: <1 initiating partons found in event " << event << endl; 
        return return_false();
    }

    IP_jet = make_pseudojet((*ip_pt)[0], (*ip_eta)[0], (*ip_phi)[0], (*ip_E)[0]);

    if (IP_jet.eta()< IP_eta_min || IP_jet.eta() > IP_eta_max) {
        return false;
    }

    auto n_const = (*pt).size();
    // vector<fastjet::PseudoJet> parts;
    all_constituents.clear();

    for (auto i = 0; i<n_const; ++i) {
        all_constituents.push_back(make_pseudojet((*pt)[i], (*eta)[i], (*phi)[i], (*E)[i]));
        all_constituents.back().set_user_index(-1000);
    }
    auto jets = clusterer(all_constituents);
    if (jets.size() == 0) return return_false();

    // do jet matching
    select_circle.set_reference(IP_jet);
    auto closest_jets = select_circle(jets);
    has_matched_IP = (closest_jets.size()>0);
    if (!has_matched_IP) {
        dR=-1.;
        return true;
    }

    // take the highest pT jet in this distance:
    int i_match = 0;
    float pt_match = closest_jets[0].pt();
    for (int i = 1; i<closest_jets.size(); ++i) {
        if (closest_jets[i].pt() > pt_match) {
            i_match = i;
            pt_match = closest_jets[i].pt();
        }
    }
    dR = IP_jet.delta_R(closest_jets[i_match]);
    truth_jet = closest_jets[i_match];
    truth_jet_constituents = fastjet::sorted_by_pt(truth_jet.constituents());
    int index = 0;
    for (auto& C : truth_jet_constituents) {
        C.set_user_index(index++);
    }
    return true;
}
