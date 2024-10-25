#include "BkgGen.h"

using namespace fastjet;
using namespace std;

void BkgGen::init(TTree* _tree, bool _separate_bulk) 
{
    tree = _tree;
    separate_bulk = _separate_bulk;
  is_init = true;
  if (separate_bulk) {
    // read the bulk particles from a hydro
    num_bulk_events = tree->GetEntries();
  }
  tree->SetBranchStatus("bulk_*", 1);

  TBranch* test_bulk = tree->GetBranch("bulk_PID");
  if (test_bulk == nullptr) {
    std::cout << " Failing to run file because there is no file for the bulk particles" << std::endl;
    assert(false);
  }

  // if (read_background_input) {
  tree->SetBranchAddress("bulk_PID", &bulk_PID);
  tree->SetBranchAddress("bulk_pt",  &bulk_pt);
  tree->SetBranchAddress("bulk_eta", &bulk_eta);
  tree->SetBranchAddress("bulk_phi", &bulk_phi);
  tree->SetBranchAddress("bulk_E",   &bulk_E);
  
  // } else {
  // rng.SetSeed(seed);
  // fpt = new TF1("fpt",Form("x*exp(-x/%10.8f)",T),minPtCut,10.);
  // double pt_ratio = fpt->Integral(minPtCut, 10.) / fpt->Integral(0., 10.);
  // nParticles = static_cast<int> ( 
      // 2 * maxEta * dNdEta * (include_neutral ? 1.5 : 1.) * pt_ratio );
  // } 
}

vector<PseudoJet> BkgGen::operator()()
{
  assert(is_init);
  std::vector<fastjet::PseudoJet> eventBkg;

  if (separate_bulk) {
    if (i_bulk_event == num_bulk_events) {
      i_bulk_event = 0;
    }
    tree->GetEntry(i_bulk_event);
    ++i_bulk_event;
  }

  int n_const = bulk_PID->size();
  for (int i = 0; i < n_const; ++i)
  {
    fastjet::PseudoJet part;
    const auto pt = (*bulk_pt)[i];
    const auto eta = (*bulk_eta)[i];
    const auto phi = (*bulk_phi)[i];

    part.reset_momentum(pt * cos(phi), pt * sin(phi), pt * sinh(eta), (*bulk_E)[i]);
    /* part.reset_PtYPhiM((*bulk_pt)[i],(*bulk_eta)[i],(*bulk_phi)[i],(*bulk_E)[i]); */
    part.set_user_index(-1); // -1 are background particles
    eventBkg.push_back(part);
  }
  // } else { 
  //     for (int i=0;i<nParticles;i++)
  //     {
  //       const double eta { rng.Uniform(-maxEta,maxEta) };
  //       const double phi { rng.Uniform(-M_PI,M_PI)     };
  //       const double pt  { fpt->GetRandom()            };

  //       fastjet::PseudoJet part;
  //       part.reset_PtYPhiM(pt,eta,phi,PionMass);
  //       part.set_user_index(-1); // -1 are background particles

  //       /* part.set_user_index(  // set user_index for charge: -1, 0, or 1 */
  //           /* (rng.Uniform(0,1) > chargedRatio) ? 100 */
  //           /* : rng.Integer(2) == 0 ? -101 : 101 ); */

  //       eventBkg.push_back(part);
  //     }
  // }
  return eventBkg;
}
