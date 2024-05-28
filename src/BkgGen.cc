#include "BkgGen.h"

using namespace fastjet;
using namespace std;

void BkgGen::init() {
  is_init = true;
  rng.SetSeed(seed);
  fpt = new TF1("fpt",Form("x*exp(-x/%10.8f)",T),minPtCut,10.);
  double pt_ratio = fpt->Integral(minPtCut, 10.) / fpt->Integral(0., 10.);
  nParticles = static_cast<int> ( 
      2 * maxEta * dNdEta * (include_neutral ? 1.5 : 1.) * pt_ratio );
}

vector<PseudoJet> BkgGen::operator()()
{
  assert(is_init);
  std::vector<fastjet::PseudoJet> eventBkg;

  for (int i=0;i<nParticles;i++)
  {
    const double eta { rng.Uniform(-maxEta,maxEta) };
    const double phi { rng.Uniform(-M_PI,M_PI)     };
    const double pt  { fpt->GetRandom()            };

    fastjet::PseudoJet part;
    part.reset_PtYPhiM(pt,eta,phi,PionMass);
    part.set_user_index(-1); // -1 are background particles

    /* part.set_user_index(  // set user_index for charge: -1, 0, or 1 */
        /* (rng.Uniform(0,1) > chargedRatio) ? 100 */
        /* : rng.Integer(2) == 0 ? -101 : 101 ); */

    eventBkg.push_back(part);
  }
  return eventBkg;
}
