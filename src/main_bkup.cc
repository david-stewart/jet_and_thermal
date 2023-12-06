/*
   This file will do the following:

   To do:
        edit pTHatMin range for LHC to run from 0 to 7TeV (just like for RHIC runs [0,200] GeV
   Edit the input parameters:
   input:
*/


/* #include "EA_set.h" */
#include "TMath.h"
#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <array>

#include "TRandom3.h"
#include <sstream>
#include <iostream>
#include <cstring>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
/* #include "Pythia8Plugins/FastJet3.h" */
#include "Math/ProbFunc.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;


int main() {
    TFile* f_out = new TFile ("main.root","recreate");
    const double sNN { 7000 };
    const double pTHatMin { 1000. };
    const double pTHatMax { 2000. };

    string name_type = "pp";

    int idA; int idB;
    if      (name_type == "pp" )  { idA = 2212; idB = 2212; }
    else if (name_type == "pAu")  { idA = 2212; idB = 1000822080; }
    else if (name_type == "AuAu") { idA = 1000822080; idB = 1000822080; }

    TRandom3 rand{0};

    // Select FastJet parameters
    const double R_FULL       = 0.4;    // Jet size.
    const double R_CH         = 0.4;    // Jet size.
    const double pTMin        = 0.2;    // Min jet pT.
    Pythia pythia;

    const double etaMax { 10. };

    const int    nbins_pt { 400 };
    const double lobin_pt { 0. };
    const double hibin_pt {100. * sNN/200 };
    TH1D* kshort { new TH1D("kshort","K0 short multiplicity in |eta|<0.5;#it{p}_{T};dK/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* kaon   { new TH1D("kaon",  "kaon multiplicity in |eta|<0.5;#it{p}_{T};dK^{-}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* antikaon   { new TH1D("antikaon",  
            "antikaon multiplicity in |eta|<0.5;#it{p}_{T};dK/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* pion   { new TH1D("pion",  "pion multiplicity in |eta|<0.5;#it{p}_{T};d#pi/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* antipion   { new TH1D("antipion",  
            "anti-pion multiplicity in |eta|<0.5;#it{p}_{T};d#pi^{-}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* proton   { new TH1D("proton",  "proton multiplicity in |eta|<0.5;#it{p}_{T};d#p/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* pbar   { new TH1D("pbar", 
            "#bar{p} multiplicity in |eta|<0.5;#it{p}_{T};d#bar{p}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };

    
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
            Form("Random:seed = %i",rand.Integer(10000000)),
        }) pythia.readString(str.c_str());

    pythia.init();
    /* Pythia8::Info& info = pythia.info ; */

    fastjet::JetDefinition jetdef (fastjet::antikt_algorithm, R_FULL);
    fastjet::Selector jetrap  = fastjet::SelectorAbsRapMax(etaMax)  ; // for full jets

    cout << " Starting the pythia loop. " << endl;

    int n_events{1000};
    const double mips_min_p { 0.2 };
    for (int iEvent{0}; iEvent < n_events; ++iEvent) {
        if (!(iEvent % 10000)) {
            cout << "Finished event " << iEvent << endl;
        }
        if (!pythia.next()) continue;
        cout << " ok " << endl;

        Event& event = pythia.event;
        std::vector <fastjet::PseudoJet> part_FULL;
        std::vector <fastjet::PseudoJet> part_CH;
        float max_bemc_Et = 0.;

        for (int i {0}; i < event.size(); ++i) {
            if (!event[i].isFinal()) continue;
            auto& e { event[i] }; 
            if (!e.isVisible()) continue;

            double pAbs   { e.pAbs() };
            if (pAbs < mips_min_p) continue;

            double eta { e.eta() };
            double pt  { e.pT()  };
            /* PseudoJet p {  e.px(), e.py(), e.pz(), e.pAbs() }; // very lazy way to get psuedorapidity */
            /* float pt { p.perp() }; */
            /* float eta { static_cast<float>(p.eta()) }; // not same as rapidity because I have set masses to zero in the pseudojet */
            if ( fabs(eta) > etaMax ) continue;
            if ( fabs(eta) < 0.5 ) {
                int id { e.id() };
                if      (id ==  310 ) kshort  ->Fill(pt);
                else if (id ==  321 ) kaon    ->Fill(pt);
                else if (id == -321 ) antikaon->Fill(pt);
                else if (id ==  2212) proton  ->Fill(pt);
                else if (id == -2212) pbar    ->Fill(pt);
                else if (id ==  211 ) pion    ->Fill(pt);
                else if (id == -211 ) antipion->Fill(pt);
            }

            PseudoJet p {  e.px(), e.py(), e.pz(), e.pAbs() }; // very lazy way to get psuedorapidity
            part_FULL.push_back(p);
            if (e.isCharged()) {
                part_CH.push_back(p);
            } else {
                double et { pt / TMath::CosH(eta) };
                if (fabs(eta)<1 && et > max_bemc_Et) { max_bemc_Et = et; };
            }
        }
        
        // make and use full jets
        fastjet::ClusterSequence clustSeq_FULL(part_FULL, jetdef);
        vector<fastjet::PseudoJet> jetsFULL = sorted_by_pt(jetrap(clustSeq_FULL.inclusive_jets(pTMin)));
        // do somethign with these jets...

        // make and use charged jets
        fastjet::ClusterSequence clustSeq_CH(part_CH, jetdef);
        vector<fastjet::PseudoJet> jetsCH = sorted_by_pt(jetrap(clustSeq_CH.inclusive_jets(pTMin)));
        // do somethign with these jets...
    }
    cout << " finished PYTHIA loop " << endl;

    f_out->Write();
    f_out->Save();
    f_out->Close();

    return 0;
};
