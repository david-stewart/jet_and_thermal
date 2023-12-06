#include "RhoMedianBkgClusterer.h"

using namespace fastjet;
using namespace std;

void RhoMedianBkgClusterer::init() {
  is_init = true;
  jet_def_bkgd = fastjet::JetDefinition (algorithm, jet_R, reco_scheme, fastjet::Best); // <--
  jet_selector = !fastjet::SelectorNHardest(njets_remove) * ( fastjet::SelectorAbsEtaMax(jet_max_rap) && !fastjet::SelectorIsPureGhost() ); 
  area_def = fastjet::AreaDefinition(
      fastjet::active_area_explicit_ghosts, 
      fastjet::GhostedAreaSpec(ghost_max_rap, 1, ghost_R));
  bk_area_est = fastjet::JetMedianBackgroundEstimator(jet_selector, jet_def_bkgd, area_def);
}

float RhoMedianBkgClusterer::operator()(vector<PseudoJet>& parts) {
  assert (is_init);
  bk_area_est.set_particles(parts);
  return bk_area_est.rho();
}
