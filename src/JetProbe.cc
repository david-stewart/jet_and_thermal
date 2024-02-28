#include "JetProbe.h"
#include <iostream>
#include "TRandom3.h"
#include "Pythia8/Info.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

vector<PseudoJet> JetProbe::operator()() {
    eta = rnd3.Uniform(eta_range.first, eta_range.second);
    phi = rnd3.Uniform(0, 2*M_PI);
    jet.reset_PtYPhiM(pt, eta, phi);
    return {jet};
};
