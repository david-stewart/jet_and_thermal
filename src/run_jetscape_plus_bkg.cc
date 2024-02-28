/*
   This file will do the following:

  Do:
  + Generate truth jet (T) (PYTHIA8) particles
    + Cluster into Truth jets
  + Generate "bkg" particles (B)
    + Add to truth partcles
    + Cluster into "Measured" jets
  + Match truth and embedded jets

  // same as src/main, but reading P8particles from an input

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

#include <map>
#include <string>

using std::map;
using std::string;

// Make truth and reco jets, and put matched jets into the TTree
// also match to a qg, and see if they matched

using namespace Pythia8;
using namespace fastjet;
using namespace std;

using std::cout;
using std::endl;

int main(int nargs, char** argv) {
  /*
   * arguments:
   *   1: number of events
   *   2: ofile name
   *   3: pythia seed
   *   4: background seed
   */

    // run pthat bins:
    // 10 14
    // 14 18
    // 18 22
    // 22 24
    // 24 28
    // 28 30 
    // 30 34 
    // 34 36 
    // 36 40

  // running parmaeters
  /* int n_events        { (nargs>1) ? atoi(argv[1]) : 1      }; */
  string f_in         { (nargs>1) ? argv[1] : "test_maketupJsHadrons.root" };
  string f_outstem    { (nargs>2) ? argv[2] : "test_p8read" };
  /* double pthatmin     { (nargs>3) ? atof(argv[3]) : 0.  }; */
  /* double pthatmax     { (nargs>4) ? atof(argv[4]) : 0.  }; */
  int scropt          { (nargs>5) ? atoi(argv[5]) : 0  };
  bool print_scramble { (nargs>6) ? static_cast<bool>(atoi(argv[6])) : false  };
  int bkg_seed        { (nargs>7) ? atoi(argv[7]) :  -100  }; // defaults to equal pythia seed
  int pyth_seed       { (nargs>8) ? atoi(argv[8]) :  -100  }; // defaults to equal pythia seed

  if (pyth_seed == -100) pyth_seed = 0;
  if (bkg_seed == -100) bkg_seed = pyth_seed;

  /* ofstream fout_t::ScrOpt::Remove };xt; */
  /* fout_txt.open((f_outstem+".txt").c_str()); */

  const float MIN_TRUTH_LEAD_JET = 8.;
  bool add_background  = true;
  bool print_CA           { false };
  bool print_Pyth_and_Bkg { false };
  bool print_matched      { false };
  bool print_Pythia8_jets { false };

  TFile* fout = new TFile ((f_outstem).c_str(),"recreate");
  // the output TTree
  TTree* tree = new TTree("T", "Truth-Reco (+gluon) matches");
  
  // quarks/gluons -- only filled when matched to jets
  float qg_pt    {};
  float qg_eta   {};
  float qg_phi   {};
  int   qg_id    {};
  float rho_bkg_est    {};
  float Xsec {};
  float XsecSigma {};

  float rho_med {};

  /* tree->Branch("qg_pt",      &qg_pt); */
  /* tree->Branch("qg_eta",     &qg_eta); */
  /* tree->Branch("qg_phi",     &qg_phi); */
  /* tree->Branch("qg_id",     &qg_id); */
    
  tree->Branch("Xsec",        &Xsec);
  tree->Branch("XsecSigma",   &XsecSigma);
  tree->Branch("rho_bkg_est", &rho_bkg_est);

  // Truth jets
  JetBranch jb_pyth {};
  jb_pyth.add_to_ttree(tree, "pyth_");

  // Reco jets
  JetBranch jb_reco {{ TAG_BKGD, AREA, ANGULARITY }};
  jb_reco.nbranch_const = 10;
  jb_reco.add_to_ttree(tree, "reco_");

  TRandom3 _rand { static_cast<unsigned int>(bkg_seed) };
  PtScrambler scrambler { jb_reco.vec_lead_const_pt, _rand,
      static_cast<PtScrambler::ScrOpt>(scropt), print_scramble };


  // Background particles
  BkgGen bkgmaker;
  bkgmaker.seed            = (unsigned int) bkg_seed;
  bkgmaker.T               = 0.291;
  bkgmaker.include_neutral = true;
  bkgmaker.minPtCut        = 0.2;
  bkgmaker.init();

  // Pythia 8 particles
  /* P8Gen p8maker {}; */
  /* p8maker.seed            = pyth_seed; */
  /* p8maker.name_type       = "pp"; */
  /* p8maker.sNN             = 200.; */
  /* p8maker.pTHatMin        = pthatmin; */
  /* p8maker.pTHatMax        = pthatmax; */
  /* p8maker.collect_neutral = true; */
  /* p8maker.collect_charged = true; */
  /* p8maker.init(); */

  P8TupReader p8reader(f_in);

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

  bool last_event = false;
  int nev = 0;
  while (last_event == false) {
      if (nev % 1000 == 0) cout << " Finished " << nev << " events" << endl;
      nev++;
    // PYTHIA8 jets
    vector<PseudoJet> part_P;
    vector<PseudoJet> jets_P;

    int n_attempt = 0;
    while(true)  {
        last_event = !(p8reader.next());
        if (last_event) break;
        auto _part_P = p8reader(); // vector<fastjet::PseudoJet>
        auto _jets_P = clusterer(_part_P); // reconstructed ("truth") jets
        if ((_jets_P.size()>0) && (_jets_P[0].perp() > MIN_TRUTH_LEAD_JET)) {
            part_P = std::move(_part_P);
            jets_P = std::move(_jets_P);
            break;
        } 
        /* else { */
            /* if (_jets_P.size() == 0) std::cout << " No leading truth jet" << std::endl; */
            /* else std::cout << " Small truth jet: " << _jets_P[0].perp() << std::endl; */
        /* } */
        if (n_attempt++ > 100) {
            std::cout << " failed finding lead jet in 10000 tries. Terminating programm" << std::endl;
            return 0;
        }
    }
    if (last_event) break;

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
    auto part_M = bkgmaker(); // vector<fastjet::PseudoJet>
    part_M.insert(part_M.end(), part_P.begin(), part_P.end());

    if (false) {
        int i =0;
        cout << endl;
        for (auto& jet : part_M) {
            cout << Form("reco: %3i: pt:eta:phi ( %5.2f, %5.2f, %5.2f)", i++, jet.perp(), jet.eta(), jet.phi()) << endl;
        }
        cout << endl;
    }
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

    i_matcher(jets_P, jets_M);
    // Match the jets_P with the qg jets (quark gluon) jets
    /* i_matcher(p8maker.e5and6, jets_P); */
    bool is_first = true;
    for (auto im : i_matcher.i_matched) {
        if (is_first) {
            rho_bkg_est = bkg_est(part_M);
        }
        jb_pyth.fill(jets_P[im.first], rho_bkg_est, is_first);
        jb_reco.fill(jets_M[im.second],rho_bkg_est, is_first);
        if (is_first) is_first = false;
        // if match to pythia truth, add the gluon
        /* auto qg = p8maker.e5and6; */
        /* qg_pt  =-100.; */
        /* qg_eta =-100.; */
        /* qg_phi =-100.; */
        /* qg_id  =-1000; */
        /* bool match_first  = (qg.size()>0 && (qg[0].squared_distance(jets_P[im.second]) <=0.16)); */
        /* bool match_second = (qg.size()>0 && (qg[1].squared_distance(jets_P[im.second]) <=0.16)); */
        /* if (match_first) { */
        /*     if (match_second) { */
        /*     cout << " WARNING: truth jet matches both initializing partons! " << endl */
        /*          << "    Therefore using only the first parton." << endl; */
        /*     } */
        /*     qg_pt = qg[0].perp(); */
        /*     qg_phi = qg[0].phi(); */
        /*     qg_eta = qg[0].eta(); */
        /*     qg_id  = qg[0].user_index(); */
        /* } else if (match_second) { */
        /*     qg_pt = qg[1].perp(); */
        /*     qg_phi = qg[1].phi(); */
        /*     qg_eta = qg[1].eta(); */
        /*     qg_id  = qg[1].user_index(); */
        /* } */
        Xsec = *p8reader.Xsec;
        XsecSigma = *p8reader.XsecSigma;

        scrambler();
        tree->Fill();
    }

    /* if (print_matched) { */ 
    /*   i_matcher(jets_P, jets_M); */
    /*   cout << " --- matched jets --- " << endl; */
    /*   for (auto im : i_matcher.i_matched) { */
    /*     int P = im.first; */
    /*     int M = im.second; */
    /*     auto& _P = jets_P[P]; */
    /*     auto& _B = jets_M[M]; */

    /*     float delta_phi = _B.phi()-_P.phi(); */
    /*     while (delta_phi > M_PI)  delta_phi -= 2*M_PI; */
    /*     while (delta_phi < -M_PI) delta_phi += 2*M_PI; */

    /*     cout << Form("PJet[%2i] (%5.2f,%5.2f,pT:%5.2f) MJet[%2i] (%5.2f,%5.2f,pT:%5.2f) Delta(%6.3f,%6.3f,pT:%6.2f,R:%6.3f(w/rap) or %6.3f(w/eta))", */
    /*       P,  _P.pseudorapidity(), _P.phi(), _P.perp(), */
    /*       M, _B.pseudorapidity(), _B.phi(), _B.perp(), */
    /*       _B.pseudorapidity() - _P.pseudorapidity(), */
    /*       delta_phi, */
    /*       _B.perp()-_P.perp(), */
    /*       _P.delta_R(_B), */ 
    /*       pow(pow(_P.pseudorapidity()-_B.pseudorapidity(),2.)+pow(delta_phi,2.),0.5) */
    /*       ) << endl; */
    /*   } */
    /* } */
  }

  // see https://pythia.org/latest-manual/CrossSectionsAndWeights.html
  map<string,double> param;
  /* param["n_events"] = n_events; */
  /* param["sigmaGen"] = p8maker.pythia.info.sigmaGen(); */
  /* param["p8_n_events"] =  p8maker.pythia.info.weightSum(); */
  /* param["Xsec"] = p8maker.pythia.info.sigmaGen() / p8maker.pythia.info.weightSum(); */
  /* param["pthatmin"] = pthatmin; */
  /* param["pthatmax"] = pthatmax; */

  fout ->cd();
  fout->WriteObject(&param,"parameters");

  fout->Write();
  fout->Save();
  fout->Close();

  return 0;
};
