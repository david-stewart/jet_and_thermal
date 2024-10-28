#include "JetTree.h"
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

JetTree::JetTree(string name, bool _use_SD, SD_criteria crit) 
  : use_SD  { _use_SD}
  , SD_crit { crit }
{
  tree = new TTree(name.c_str(), "Quark Gluon Matches");
}

void JetTree::init_branches() {
  tree->Branch("qg_pt",      &qg_pt);
  tree->Branch("qg_eta",     &qg_eta);
  tree->Branch("qg_phi",     &qg_phi);
  tree->Branch("qg_id",      &qg_id);

  tree->Branch("jet_pt",     &jet_pt);
  tree->Branch("jet_eta",    &jet_eta);
  tree->Branch("jet_phi",    &jet_phi);
  tree->Branch("jet_charge", &jet_charge);
  tree->Branch("jet_nconst", &jet_nconst);

  if (use_SD) {
    tree->Branch("jet_zg", &jet_zg);
    tree->Branch("jet_Rg", &jet_Rg);
    tree->Branch("jet_mu", &jet_mu);
  }

}

void JetTree::add_jetpair(fastjet::PseudoJet& qg_jet, fastjet::PseudoJet& reco_jet) {
  qg_pt  .push_back( qg_jet.perp()       );
  qg_eta .push_back( qg_jet.eta()        );
  qg_phi .push_back( qg_jet.phi()        );
  qg_id  .push_back( qg_jet.user_index() );

  jet_pt  .push_back( reco_jet.perp() );
  jet_eta .push_back( reco_jet.eta() );
  jet_phi .push_back( reco_jet.phi() );

  int _jet_charge = 0;
  int _jet_nconst = 1;

  if (reco_jet.has_constituents()) {
    auto constituents = reco_jet.constituents();
    _jet_nconst = (int) constituents.size();
    for (auto& C : reco_jet.constituents()) {
      _jet_charge += (int) C.user_index();
    }
  } else {
    _jet_charge = (int) reco_jet.user_index();
  }

  jet_charge .push_back( _jet_charge );
  jet_nconst .push_back( _jet_nconst );

  if (use_SD) {
    fastjet::contrib::SoftDrop sd(SD_crit.beta, SD_crit.zg);
    fastjet::PseudoJet sd_jet = sd(reco_jet);
    jet_zg.push_back( sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry());
    jet_Rg.push_back( sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R() );
    jet_mu.push_back( sd_jet.structure_of<fastjet::contrib::SoftDrop>().mu()      );
  }
}

void JetTree::clear_vectors() {
  qg_pt.clear();
  qg_eta.clear();
  qg_phi.clear();
  qg_id.clear();

  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_charge.clear();
  jet_nconst.clear();

  if (use_SD) {
    jet_zg.clear();
    jet_Rg.clear();
    jet_mu.clear();
  }
}

void JetTree::fill() {
  tree->Fill();
  clear_vectors();
}

void JetTree::write() { tree->Write(); };
