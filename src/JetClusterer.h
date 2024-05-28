#ifndef JetClusterer_h
#define JetClusterer_h

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"


struct JetClusterer {
  // take in vectors of PseudoJets and cluster into final Jets
  fastjet::JetAlgorithm algorithm { fastjet::JetAlgorithm::antikt_algorithm };
  fastjet::RecombinationScheme reco_scheme { fastjet::RecombinationScheme::E_scheme };

  // calculate areas with ghosts?
  bool  calc_area     { false };
  float ghost_max_rap { 1.8   };
  float ghost_R       { 0.01   };

  // jets
  float jet_max_rap   { 1.0  };
  float jet_R         { 0.4  };

  float min_jet_pt    { 0.2  };

  // parameters to initialize
  fastjet::JetDefinition jet_def;
  fastjet::AreaDefinition area_def;
  fastjet::Selector jet_selector;

  fastjet::ClusterSequenceArea* clus_seq_area { nullptr };
  fastjet::ClusterSequence*     clus_seq      { nullptr };

  bool is_init { false };
  int  user_index { 0 };
  void init();
  std::vector<fastjet::PseudoJet> operator() (std::vector<fastjet::PseudoJet>& p0);
};

#endif
