#ifndef RhoMedianBkgClusterer_h
#define RhoMedianBkgClusterer_h

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

struct RhoMedianBkgClusterer {
  float ghost_max_rap { 1.2 };
  float ghost_R       { 0.01 };
  int   njets_remove  { 2   };
  fastjet::JetAlgorithm algorithm { fastjet::JetAlgorithm::kt_algorithm };
  fastjet::RecombinationScheme reco_scheme { fastjet::RecombinationScheme::E_scheme };

  float jet_R         { 0.4 };
  float jet_max_rap   { 0.7   };

  fastjet::JetDefinition jet_def_bkgd;
  fastjet::Selector      jet_selector;
  fastjet::AreaDefinition area_def; 
  fastjet::JetMedianBackgroundEstimator bk_area_est;

  void init();
  bool is_init { false };
  float operator()(std::vector<fastjet::PseudoJet>& parts);
};





#endif
