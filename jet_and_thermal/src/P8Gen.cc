#include "P8Gen.h"
#include <iostream>
#include "TRandom3.h"
#include "Pythia8/Info.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;


void P8Gen::init() {
  is_init = true;

  int idA, idB;
  if      (name_type == "pp" )  { idA = 2212; idB = 2212; }
  else if (name_type == "pAu")  { idA = 2212; idB = 1000822080; }
  else if (name_type == "AuAu") { idA = 1000822080; idB = 1000822080; }

  for (auto str : vector<string>{ 
      Form("Beams:eCM = %f",sNN),  // run at 7 GeV for now
      "HardQCD:all = on",
      Form("PhaseSpace:pTHatMin = %f",pTHatMin),
      Form("PhaseSpace:pTHatMax = %f",pTHatMax),
      Form("Beams:idA = %i", idA), // moving in +z direction
      Form("Beams:idB = %i", idB), // moving in +z direction
      /* "Beams:idB = 1000822080", // moving in -z direction AuAu */
      /* "Beams:idB = 2212", // moving in -z direction AuAu */
                             // Au ion (from http://home.thep.lu.se/~torbjorn/pythia82html/BeamParameters.html#section0 )
      "Random:setSeed = on",
      "ParticleDecays:limitRadius = on",
      Form("ParticleDecays:rMax = %f", 10.),
      Form("Random:seed = %i",seed)
  }) pythia.readString(str.c_str());
  pythia.init();

  e5and6 .push_back(fastjet::PseudoJet());
  e5and6 .push_back(fastjet::PseudoJet());

};

vector<PseudoJet> P8Gen::operator()() {
  assert(is_init);

  int ntries = 1;
  bool good_event = pythia.next();
  while (!good_event) {
    ++ntries;
    if (ntries > nMaxBadGenerations) {
      cout << " Pythia8 failed to generate a good event " << ntries << " successive times." << endl;
      cout << " Aborting program." << endl;
      assert(ntries <= nMaxBadGenerations);
    }
    good_event = pythia.next();
  } 

  /* auto info = pythia.info; */
  /* weightSum = pythia.info.weightSum(); */
  /* sigmaGen  = pythia.info.sigmaGen(); */
  /* cout << " sigmaGen: " << (sigmaGen*weightSum) << " and weightSum: " << weightSum << Form(" -> %10.5g",sigmaGen) << endl; */

  auto& event = pythia.event ;

  std::vector <fastjet::PseudoJet> p_vec; // charged vec

  if (event.size() < 7) return {};

  // get e5 and e6 (the scattered particles)
  auto& e = event[5];
  e5and6[0].reset_PtYPhiM(e.pT(), e.eta(), e.phi(), e.m());
  e5and6[0].set_user_index(e.id());

  e = event[6];
  e5and6[1].reset_PtYPhiM(e.pT(), e.eta(), e.phi(), e.m());
  e5and6[1].set_user_index(e.id());

  for (int i {0}; i < event.size(); ++i) {
    auto& e     { event[i] };
    if (!event[i].isFinal()) continue;

    bool isCharged = e.isCharged();
    if ((isCharged  && !collect_charged) 
    ||  (!isCharged && !collect_neutral)) continue;

    double pAbs { e.pAbs() };

    double eta { e.eta() };
    if ( fabs(eta) > maxEta ) continue;

    double pt { e.pT() };
    if (pt < minPtCut) continue;

    double phi { e.phi() };

    fastjet::PseudoJet p;
    p.reset_PtYPhiM(pt, eta, phi, usePionMass ? PionMass : e.m());
    p.set_user_index(static_cast<int>(e.charge()));
    p_vec.push_back(p);
  }
  return p_vec;
};

P8TupReader::P8TupReader(const std::string& fname) 
    : fin { new TFile(fname.c_str(), "read") }
    , event_reader { "events", fin }
    , part_reader  { "particles", fin }
    , px    { part_reader, "px" }
    , py    { part_reader, "py" }
    , pz    { part_reader, "pz" }
    , E     { part_reader, "E" }
    , nPart { event_reader, "nhadrons" }
    , Xsec  { event_reader, "Xsec" }
    , XsecSigma { event_reader, "XsecSigma" }
{ };

bool P8TupReader::next() {
    return event_reader.Next();
}

std::vector<fastjet::PseudoJet> P8TupReader::operator()() {
    int npart = static_cast<int>(*nPart);
    std::vector <fastjet::PseudoJet> p_vec; // charged vec
    for (int i=0;i<npart;++i) {
        bool has_next = part_reader.Next();
        if (!has_next) {
            std::cout << "Fatal error -- trying to read too many input hadrons!" << std::endl;
            break;
        }
        fastjet::PseudoJet p { *px, *py, *pz, *E };
        p.set_user_index(i); // jetscape particles are 2
        p_vec.push_back(p);
        p.set_user_index(static_cast<int>(-100));
    }
    return p_vec;
}

std::tuple<bool, fastjet::PseudoJet, std::vector<fastjet::PseudoJet>> 
P8TupReader::lead_jet_only(JetClusterer& clusterer) 
{
    auto part = this->operator()();
    auto jets = clusterer(part);
    if (jets.size()==0) { 
        return std::make_tuple(false, fastjet::PseudoJet(), vector<fastjet::PseudoJet>());
    } else {
      auto parts = fastjet::sorted_by_pt(jets[0].constituents());
      int n=0;
      for (auto& C : parts) {
        C.set_user_index(n++);
      }
      return std::make_tuple(true, jets[0], parts);
    }
}

