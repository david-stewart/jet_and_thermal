#include "JetClusterer.h"

using namespace fastjet;
using namespace std;

#include <cassert>

void JetClusterer::init() {
  is_init = true;
  jet_def = fastjet::JetDefinition(algorithm, jet_R, reco_scheme, fastjet::Best);
  fastjet::GhostedAreaSpec area_spec(ghost_max_rap);
  area_def = (active_area, area_spec);
  /* if (calc_area) area_def = fastjet::AreaDefinition( */
      /* fastjet::active_area_explicit_ghosts, */ 
      /* fastjet::GhostedAreaSpec(ghost_max_rap, 1, ghost_R)); */
  jet_selector = fastjet::SelectorAbsEtaMax(jet_max_rap) && !fastjet::SelectorIsPureGhost(); 
}

vector<PseudoJet> JetClusterer::operator()(vector<PseudoJet>& pseudojets)
{
  assert (is_init);
  if (calc_area) {
    if (clus_seq_area != nullptr) delete clus_seq_area;
    clus_seq_area = new fastjet::ClusterSequenceArea(pseudojets, jet_def, area_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt( jet_selector( clus_seq_area->inclusive_jets(min_jet_pt) ));
    return jets;
  } else {
    if (clus_seq != nullptr) delete clus_seq;
    clus_seq = new fastjet::ClusterSequence (pseudojets, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt( jet_selector( clus_seq->inclusive_jets(min_jet_pt) ));
    return jets;
  }
}
