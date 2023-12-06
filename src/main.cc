/*
   This file will do the following:

  Do:
  + Generate truth jet (T) (PYTHIA8) particles
    + Cluster into Truth jets
  + Generate "bkg" particles (B)
    + Add to truth partcles
    + Cluster into "Measured" jets
  + Match truth and embedded jets

*/

#include "BkgGen.h"
#include "P8Gen.h"
#include "JetClusterer.h"
#include "JetIndicesMatcher.h"
#include "RhoMedianBkgClusterer.h"
#include "JTWalker.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVectorD.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "JetTree.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;


int main(int nargs, char** argv) {
  /*
   * arguments:
   *   1: number of events
   *   2: ofile name
   *   3: pythia seed
   *   4: background seed
   */

  // running parmaeters
  int n_events     { (nargs>1) ? atoi(argv[1]) : 1     };
  string f_outstem { (nargs>2) ? argv[2] : "pithy"     };
  int pyth_seed    { (nargs>3) ? atoi(argv[3]) : 1     }; // -1 is default for Pythia
  int bkg_seed     { (nargs>4) ? atoi(argv[4]) :  -100 }; // defaults to equal pythia seed
  if (bkg_seed == -100) bkg_seed = pyth_seed;

  ofstream fout_txt;
  fout_txt.open((f_outstem+".txt").c_str());

  TFile* fout = new TFile ((f_outstem+".root").c_str(),"recreate");

  JetTree tree{"T",true};
  tree.init_branches();

  // print arguments:
  bool add_background    = false;
  bool print_Pythia8_jets { false };
  bool print_CA { true };
  bool print_Pyth_and_Bkg { false };
  bool print_matched      { false };

  
  // Background particles
  BkgGen bkgmaker;
  bkgmaker.seed            = (unsigned int) bkg_seed;
  bkgmaker.T               = 0.291;
  bkgmaker.include_neutral = true;
  bkgmaker.minPtCut        = 0.2;
  bkgmaker.init();

  // Pythia 8 particles
  P8Gen p8maker {};
  p8maker.seed            = pyth_seed;
  p8maker.name_type       = "pp";
  p8maker.sNN             = 200.;
  p8maker.pTHatMin        = 30.;
  p8maker.pTHatMax        = 80.;
  p8maker.collect_neutral = true;
  p8maker.collect_charged = true;
  p8maker.init();

  // JetClusterer
  JetClusterer clusterer {};
  clusterer.calc_area = false;
  clusterer.min_jet_pt = 5.;
  clusterer.init();


  // JetCA_clusterer
  JetClusterer clusterer_CA {};
  clusterer_CA.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  clusterer_CA.jet_max_rap = 2.; // shouldn't be necessary
  clusterer_CA.jet_R       = 2.; // shouldn't be necessary
  clusterer_CA.init();

  // i_matcher for jets
  JetIndicesMatcher i_matcher { 0.4 }; 

  // get background rho estimation
  RhoMedianBkgClusterer bkg_est{};
  bkg_est.init();

  for (int nev=0;nev<n_events;++nev) {
    // PYTHIA8 jets
    auto part_P = p8maker(); // vector<fastjet::PseudoJet>
    /* cout << " FIXME(" << nev <<") gen: " << p8maker.weightSum << " and " << p8maker.sigmaGen << endl; */
    /* cout << " FIXME weight: " << p8maker.pythia.info.weight() << endl; */
    auto jets_P = clusterer(part_P); // reconstructed ("truth") jets

    if (print_Pythia8_jets) { 
      cout << " --- Jets from PYTHIA8 --- size: " << jets_P.size() << endl;
      int i {0};
      for (auto& jet : jets_P) {
        auto comps = jet.constituents();
        cout << Form("jet[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f) nconst(%3i)",
          i++, jet.eta(),jet.phi(),jet.perp(), (int) comps.size()) << endl;
        int ic {0};
        for (auto& C : comps) {
          cout << Form("  const[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f)",
            ic++, C.eta(),C.phi(),C.perp()) << endl;

          // print out the path
        }
      }
    }

    // Add in background
    vector<fastjet::PseudoJet>& jets_M = jets_P;
    if (add_background) {
      auto part_B = bkgmaker(); // vector<fastjet::PseudoJet>
      part_B.insert(part_B.end(), part_P.begin(), part_P.end());
      jets_M = clusterer(part_B);
    }

    if (print_CA) for (auto& jet : jets_M) {
      auto comps = jet.constituents();
      auto jets_CA = clusterer_CA(comps);

      if (jets_CA.size() == 1) {
        SD_criteria criteria;
        criteria.zg = 0.1;
        criteria.min_pt = 1.0;
        criteria.beta = 0.;
        criteria.R0 = 0.;
        TreeSplitPrinter print1 { jets_CA[0], criteria, "", fout_txt };
        /* TreeSplitPrinter print2 { jets_CA[0], criteria, "", std::cout }; */
      }
      break;
    }

    if (print_Pyth_and_Bkg) { 
      cout << " --- Pythia8 Jets ";
      if (add_background) cout << " + background";
      cout << "--- " << endl;

      int i {0};
      for (auto& jet : jets_M) {
        auto comps = jet.constituents();
        auto jets_CA = clusterer_CA(comps);
        cout << Form("jet[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f) nconst(%3i) nCArecluster(%3i)",
          i++, jet.eta(),jet.phi(),jet.perp(), (int) comps.size(), (int) jets_CA.size()) << endl;
        int ic {0};
        for (auto& C : comps) {
          cout << Form("  const[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f)",
            ic++, C.eta(),C.phi(),C.perp()) << endl;
        }

        // Now recluster the jet and show what it is
        int ica {0};
        for (auto cjet : jets_CA) {
          auto ca_comps = jet.constituents();
          cout << Form("  -> jet CA[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f) nconst(%3i)",
            ica++, cjet.eta(),cjet.phi(),cjet.perp(), (int) ca_comps.size()) << endl;
          int ic {0};
          for (auto& C : ca_comps) {
            cout << Form("    -> const[%2i] eta:phi:pt(%5.2f,%5.2f,%5.2f)",
              ic++, C.eta(),C.phi(),C.perp()) << endl;
          }
        }
        
      }
    }


    // Match the jets_P with the qg jets (quark gluon) jets
    i_matcher(p8maker.e5and6, jets_P);

    for (auto im : i_matcher.i_matched) {
      int qg = im.first;
      int reco = im.second;
      tree.add_jetpair(p8maker.e5and6[im.first], jets_P[im.second]);
    }
    tree.fill();

    // Match the jets
    if (print_matched) { 
      i_matcher(jets_P, jets_M);
      cout << " --- matched jets --- " << endl;
      for (auto im : i_matcher.i_matched) {
        int P = im.first;
        int M = im.second;
        auto& _P = jets_P[P];
        auto& _B = jets_M[M];

        float delta_phi = _B.phi()-_P.phi();
        while (delta_phi > M_PI)  delta_phi -= 2*M_PI;
        while (delta_phi < -M_PI) delta_phi += 2*M_PI;

        cout << Form("PJet[%2i] (%5.2f,%5.2f,pT:%5.2f) MJet[%2i] (%5.2f,%5.2f,pT:%5.2f) Delta(%6.3f,%6.3f,pT:%6.2f,R:%6.3f(w/rap) or %6.3f(w/eta))",
          P,  _P.pseudorapidity(), _P.phi(), _P.perp(),
          M, _B.pseudorapidity(), _B.phi(), _B.perp(),
          _B.pseudorapidity() - _P.pseudorapidity(),
          delta_phi,
          _B.perp()-_P.perp(),
          _P.delta_R(_B), 
          pow(pow(_P.pseudorapidity()-_B.pseudorapidity(),2.)+pow(delta_phi,2.),0.5)
          ) << endl;
      }
    }
  }

  fout->Write();
  fout->Save();
  fout->Close();

  return 0;
};
