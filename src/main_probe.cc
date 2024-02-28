/*
   This file will do the following:

  Do:
  + Generate truth jet (T) (PYTHIA8) particles
    + Cluster into Truth jets
        scrambler();
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
/* #include "JTWalker.h" */

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVectorD.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "JetBranch.h"
#include "PtScrambler.h"
#include "JetProbe.h"


// Make truth and reco jets, and put matched jets into the TTree
// also match to a qg, and see if they matched

using namespace fastjet;
using namespace std;

using std::cout;
using std::endl;

int main(int nargs, char** argv) {
  /*
   * arguments:
   *   1: number of events
   *   2: ofile name
   *   4: background seed
   */

  // running parmaeters
  int n_events        { (nargs>1) ? atoi(argv[1]) : 1      };
  string f_outstem    { (nargs>2) ? argv[2] : "matchedjet" };
  double pthatmin     { (nargs>3) ? atof(argv[3]) : 10.  };
  double pthatmax     { (nargs>4) ? atof(argv[4]) : 20.  };
  int scropt          { (nargs>3) ? atoi(argv[3]) : 0  };
  bool print_scramble { (nargs>4) ? static_cast<bool>(atoi(argv[4])) : false  };
  int bkg_seed        { (nargs>5) ? atoi(argv[5]) :  0  }; // defaults to equal pythia seed


  /* ofstream fout_t::ScrOpt::Remove };xt; */
  /* fout_txt.open((f_outstem+".txt").c_str()); */

  const float MIN_TRUTH_LEAD_JET = 10.;
  bool add_background  = true;
  bool print_CA           { false };
  bool print_probe_and_bkg { false };
  bool print_matched      { false };


  TFile* fout = new TFile ((f_outstem+".root").c_str(),"recreate");
  // the output TTree
  TTree* tree = new TTree("T", "Truth-Reco (+gluon) matches");
  // quarks/gluons -- only filled when matched to jets
  float rho_bkg_est    {};
  float rho_med {};

  // JetProbe
  TRandom3 _rand { static_cast<unsigned int>(bkg_seed) };
  JetProbe probe {_rand, 0.7};

  JetBranch jb_probe {};
  jb_probe.add_to_ttree(tree, "probe_");

  // Reco jets
  JetBranch jb_reco {{ TAG_BKGD, AREA, ANGULARITY }};
  jb_reco.nbranch_const = 10;
  jb_reco.add_to_ttree(tree, "reco_");

  PtScrambler scrambler { jb_reco.vec_lead_const_pt, _rand,
      static_cast<PtScrambler::ScrOpt>(scropt), print_scramble };


  // Background particles
  BkgGen bkgmaker;
  bkgmaker.seed            = (unsigned int) bkg_seed;
  bkgmaker.T               = 0.291;
  bkgmaker.include_neutral = true;
  bkgmaker.minPtCut        = 0.2;
  bkgmaker.init();

  // JetClusterer
  JetClusterer clusterer {};
  clusterer.calc_area = false;
  clusterer.min_jet_pt = 5.;
  clusterer.init();

  // Jetclusterer with area
  JetClusterer clusterer_area {};
  clusterer_area.calc_area = jb_reco.fill_area;
  clusterer_area.min_jet_pt = 5.;
  clusterer_area.init();

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
      if (nev % 1000 == 0) cout << " Finished " << nev << " events" << endl;
    // probe jet
    vector<PseudoJet> part_P = probe();
    vector<PseudoJet> jets_P = clusterer(part_P);

    // Add in background
    auto part_M = bkgmaker(); // vector<fastjet::PseudoJet>
    part_M.insert(part_M.end(), part_P.begin(), part_P.end());
    auto jets_M = clusterer_area(part_M);

    if (print_CA) for (auto& jet : jets_M) {
      auto comps = jet.constituents();
      auto jets_CA = clusterer_CA(comps);

      if (jets_CA.size() == 1) {
        SD_criteria criteria;
        criteria.zg = 0.1;
        criteria.min_pt = 1.0;
        criteria.beta = 0.;
        criteria.R0 = 0.;
        /* TreeSplitPrinter print1 { jets_CA[0], criteria, "", fout_txt }; */
        /* TreeSplitPrinter print2 { jets_CA[0], criteria, "", std::cout }; */
      }
      break;
    }


    if (print_probe_and_bkg) { 
      cout << " --- JetProbe ";
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

    i_matcher(jets_M, jets_P);
    // Match the jets_P with the qg jets (quark gluon) jets
    /* i_matcher(p8maker.e5and6, jets_P); */
    bool is_first = true;
    for (auto im : i_matcher.i_matched) {
        if (is_first) {
            rho_bkg_est = bkg_est(part_M);
        }
        jb_reco.fill(jets_M[im.first],rho_bkg_est, is_first);
        jb_probe.fill(jets_P[im.second],rho_bkg_est, is_first);
        scrambler();
        tree->Fill();
        if (is_first) is_first = false;
    }

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
