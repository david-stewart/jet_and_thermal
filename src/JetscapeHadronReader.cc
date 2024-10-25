#include "JetscapeHadronReader.h"
#include <iostream>
#include "TRandom3.h"
#include "Pythia8/Info.h"
#include "TLorentzVector.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

/* std::tuple<bool, fastjet::PseudoJet, std::vector<fastjet::PseudoJet>> */ 
/* P8TupReader::lead_jet_only(JetClusterer& clusterer) */ 
/* { */
/*     auto part = this->operator()(); */
/*     auto jets = clusterer(part); */
/*     if (jets.size()==0) { */ 
/*         return std::make_tuple(false, fastjet::PseudoJet(), vector<fastjet::PseudoJet>()); */
/*     } else { */
/*       auto parts = fastjet::sorted_by_pt(jets[0].constituents()); */
/*       int n=0; */
/*       for (auto& C : parts) { */
/*         C.set_user_index(n++); */
/*       } */
/*       return std::make_tuple(true, jets[0], parts); */
/*     } */
/* } */

void JetscapeHadronReader::write() {
    LJ_deltaR ->Write();
    LJ_delta_phieta ->Write();
    Edens_deltaR ->Write();
    Edens_delta_phieta ->Write();

}

JetscapeHadronReader::JetscapeHadronReader(const std::string& fname) 
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

bool JetscapeHadronReader::next_event() {
    i_event++;
    if (i_event >= nEvents) {
        return false;
    } 
    if (print_QA) cout << " FIXME Event: " << i_event << endl;
    return tree->GetEntry(i_event);
}

fastjet::PseudoJet JetscapeHadronReader::make_pseudojet(float pt, float eta, float phi, float E) {
    TLorentzVector C;
    C.SetPtEtaPhiE(pt, eta, phi, E);
    return fastjet::PseudoJet(C);
}

std::array<double,3> JetscapeHadronReader::deltas(fastjet::PseudoJet& a, fastjet::PseudoJet& b) {
    auto delta_phi = b.phi() - a.phi();
    while (delta_phi>M_PI) delta_phi -= 2*M_PI;
    while (delta_phi<-M_PI) delta_phi += 2*M_PI;
    auto delta_eta = b.eta() - a.eta();
    auto delta_R = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
    return {delta_phi, delta_eta, delta_R};
}

bool JetscapeHadronReader::cluster_jets(JetClusterer& clusterer) {
    int n_IP = ip_pt->size();
    if (n_IP < 2) {
        if (print_QA) cout << "FIXME: <2 initiating partons found in event " << event << endl; 
        return false;
    }

    IP_leading = make_pseudojet((*ip_pt)[0], (*ip_eta)[0], (*ip_phi)[0], (*ip_E)[0]);
    IP_sub     = make_pseudojet((*ip_pt)[1], (*ip_eta)[1], (*ip_phi)[1], (*ip_E)[1]);

    auto n_jets = (*pt).size();
    vector<fastjet::PseudoJet> parts;
    for (auto i = 0; i<n_jets; ++i) { 
        parts.push_back(make_pseudojet((*pt)[i], (*eta)[i], (*phi)[i], (*E)[i]));
        std::cout << " FIXME A0: E: " << (*E)[i] << std::endl;
    }
    auto jets = clusterer(parts);
    if (jets.size() == 0) {
        has_leading_jet = false;
        leading_jet = fastjet::PseudoJet(); 
        IP_match_index = -1;
        leading_jet_constituents.clear();
        return true;
    }  
    for (auto& FIXME_jet : jets) {
        std::cout << " FIXME A1: jet-rap(" << FIXME_jet.rap() <<") eta(" << FIXME_jet.eta() <<") diff(" << (FIXME_jet.eta()-FIXME_jet.rap()) << std::endl;
    }
    leading_jet = jets[0];
    has_leading_jet = true;

    // fill in the histograms
    auto delta_lead = deltas(IP_leading, leading_jet);
    auto delta_sub =  deltas(IP_sub, leading_jet);
    auto& IP = (delta_lead[2] < delta_sub[2]) ? IP_leading : IP_sub;
    auto& delta = (delta_lead[2] < delta_sub[2]) ? delta_lead : delta_sub;
    LJ_deltaR->Fill(delta[2]);
    LJ_delta_phieta->Fill(delta[0], delta[1]);
    float area = Edens_deltaR->GetXaxis()->GetBinWidth(1) * 
                 Edens_deltaR->GetYaxis()->GetBinWidth(1);
    for (auto& C : parts) {
        auto D = this->deltas(IP, C);
        Edens_deltaR->Fill(D[2],C.pt()/(area*(1+D[2])));
        Edens_delta_phieta->Fill(D[0], D[1], C.pt());
    }

    // set the constituent indices
    leading_jet_constituents = fastjet::sorted_by_pt(leading_jet.constituents());
    int index = 0;
    for (auto& C : leading_jet_constituents) {
        C.set_user_index(index++);
    }

    // match the leading jet
    auto R0 = leading_jet.delta_R(IP_leading);
    auto R1 = leading_jet.delta_R(IP_sub);
    if (R0 < R1) {
        if (R0 < 0.3)  IP_match_index =  0;
        else           IP_match_index = -1;
    } else {
        if (R1 < 0.3)  IP_match_index =  1;
        else           IP_match_index = -1;
    }
    return true;
}
