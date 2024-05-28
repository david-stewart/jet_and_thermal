/*
   This file will do the following:

Do:
+ Get the leading initiating parton (IP)
+ Cluster the jets (R=0.4)
+ Find the closest reconstructed jets (within DeltaR<=0.3)
+ Match to the highest pT closest jet
+ Write out the leading partons to the tree
+ Fill out of THnSparse: IP_pt, jet_pT, framentation function, Xsec
*/

#include "BkgGen.h"
// #include "JetscapeHadronReader.h"
#include "JetClusterer.h"
#include "IP_geom_matcher.h"
// #include "JetIndicesMatcher.h"
#include "RhoMedianBkgClusterer.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
// #include "TVectorD.h"
// #include "fastjet/ClusterSequenceArea.hh"
// #include "fastjet/PseudoJet.hh"
// #include "fastjet/Selector.hh"
#include "JetBranch.h"
// #include "PtScrambler.h"

#include <map>
#include <string>

using std::map;
using std::string;

// Make truth and reco jets, and put matched jets into the TTree
// also match to a qg, and see if they matched

// using namespace Pythia8;
using namespace fastjet;
using namespace std;

int main(int nargs, char **argv)
{
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
  string f_in{(nargs > 1) ? argv[1] : "test_in.root"};
  string f_outstem{(nargs > 2) ? argv[2] : "test_out.root"};
  bool   fill_vec_cpt{(nargs > 3) ? static_cast<bool>(atoi(argv[3])) : false};
  const int N_BKGSAMPLES{(nargs > 4) ? atoi(argv[4]) : 1};
  int max_nevents{(nargs > 5) ? atoi(argv[5]) : -1};

  cout << " Input: " << f_in << "->" << f_outstem <<" vec("<<fill_vec_cpt<<") nemb="<<N_BKGSAMPLES<<" ev="<<max_nevents<<endl;

  /* double pthatmin     { (nargs>3) ? atof(argv[3]) : 0.  }; */
  /* double pthatmax     { (nargs>4) ? atof(argv[4]) : 0.  }; */
  // int scropt{(nargs > 5) ? atoi(argv[5]) : 0};
  // bool print_scramble{(nargs > 6) ? static_cast<bool>(atoi(argv[6])) : false};
  int bkg_seed{(nargs > 6) ? atoi(argv[6]) : -100};  // defaults to equal pythia seed
  // int pyth_seed{(nargs > 8) ? atoi(argv[8]) : -100}; // defaults to equal pythia seed

  // if (pyth_seed == -100)
    // pyth_seed = 0;
  if (bkg_seed == -100) bkg_seed = 0.;

  const float MIN_TRUTH_LEAD_JET = 0.;
  bool add_background = true;
  bool print_CA{false};
  bool print_Pyth_and_Bkg{false};
  bool print_matched{false};
  bool print_Pythia8_jets{false};

  TFile *fout = new TFile((f_outstem).c_str(), "recreate");
  TTree *tree = new TTree("T", "Trees for the initiating parton");

  IP_geom_matcher js_reader(f_in);

  float pt_IP;
  float eta_IP;
  float phi_IP;

  bool matched_TtoR  = false;  // truth to reco
  bool matched_IPtoT = false; // IP to truth

  float dR_IPtoT = -1.;
  float dR_TtoR  = -1.;

  // vector<float> pt_truth_const; // the pt of the constituents of the leading truth jet
  // vector<bool> ismatched_const; // whether the constituents of the leading truth jet are matched to a reco jet
  /* float reco_pt_truthconst;     // the collective "jet pt" of the ismatched_const constituents */

  float rho_bkg_thermal{}; // background est. without pythia
  float rho_bkg_thermalandjet{}; // background est. with added jetscape (ignore 2 leading hardest)
  float Xsec{};
  float XsecSigma{};

  float resid_rhoA{}; // using rho_bkg_thermalandjet

  tree->Branch("matched_IPtoT", &matched_IPtoT);
  tree->Branch("matched_TtoR",  &matched_TtoR);

  tree->Branch("dR_IPtoT", &dR_IPtoT);
  tree->Branch("dR_TtoR", &dR_TtoR);

  tree->Branch("IP_pt", &pt_IP);
  tree->Branch("IP_eta", &eta_IP);
  tree->Branch("IP_phi", &phi_IP);

  tree->Branch("resid_rhoA", &resid_rhoA);
  /* tree->Branch("pt_truth_const", &pt_truth_const); */
  /* tree->Branch("ismatched_const", &ismatched_const); */
  /* tree->Branch("reco_pt_truthconst", &reco_pt_truthconst); */

  tree->Branch("Xsec", &Xsec);
  tree->Branch("XsecSigma", &XsecSigma);
  tree->Branch("rho_bkg_thermal", &rho_bkg_thermal);
  tree->Branch("rho_bkg_thermalandjet", &rho_bkg_thermalandjet);

  vector<JETBRANCH_OPTIONS> truth_options = {};
  if (fill_vec_cpt) truth_options.push_back(JETBRANCH_OPTIONS::VEC_CONSTITUENTS);
  JetBranch jb_truth{truth_options};
  jb_truth.nbranch_const = 10;
  jb_truth.add_to_ttree(tree, "truth_");

  JetBranch jb_reco{{TAG_BKGD, AREA, ANGULARITY}};
  jb_reco.nbranch_const = 10;
  jb_reco.add_to_ttree(tree, "reco_");

  TRandom3 _rand{static_cast<unsigned int>(bkg_seed)};

  // Background particles
  BkgGen bkgmaker;
  bkgmaker.seed = (unsigned int)bkg_seed;
  bkgmaker.T = 0.291;
  bkgmaker.include_neutral = true;
  bkgmaker.minPtCut = 0.2;
  bkgmaker.init();

  // JetClusterer
  JetClusterer clusterer{};
  clusterer.jet_R = 0.4;
  clusterer.calc_area = false;
  clusterer.min_jet_pt = 0.;
  clusterer.init();

  // Jetclusterer with area
  JetClusterer clusterer_area{};
  clusterer_area.jet_R = 0.4;
  clusterer_area.calc_area = jb_reco.fill_area;
  clusterer_area.min_jet_pt = 0.;
  clusterer_area.init();

  // get background rho estimation
  RhoMedianBkgClusterer bkg_est{};
  bkg_est.njets_remove = 0;
  bkg_est.init();

  RhoMedianBkgClusterer bkg_est_wJetScape{};
  bkg_est_wJetScape.njets_remove = 2;
  bkg_est_wJetScape.init();

  fastjet::PseudoJet blank_constituent(0.,0.,0.,0.);
  blank_constituent.set_user_index(-1);
  vector<fastjet::PseudoJet> blank_constituents = {blank_constituent};
  auto blank_jets = clusterer(blank_constituents);

  // JetIndicesMatcher i_matcher{0.3}; // slightly tighten the matching radius
  fastjet::Selector select_circle = SelectorCircle(0.3);

  // get background rho estimation
  int nev = 0;
  while (js_reader.next_event())
  {
    if (max_nevents != -1 && js_reader.i_event > max_nevents)
      break;
    if (!js_reader.cluster_match(clusterer))
      continue; // no IP close

    if (nev % 1000 == 0)
      cout << " Finished " << nev << " events" << endl;
    nev++;

    float zero = 0.;

    Xsec = js_reader.sigma;
    XsecSigma = js_reader.sigma_error;
    pt_IP = js_reader.IP_jet.pt();
    eta_IP = js_reader.IP_jet.eta();
    phi_IP = js_reader.IP_jet.phi();

    matched_IPtoT = js_reader.has_matched_IP;
    dR_IPtoT = js_reader.dR;

    if (!matched_IPtoT)
    { 
        matched_IPtoT = false;
        matched_TtoR  = false;
        dR_TtoR = -1.;
        rho_bkg_thermal = 0.; 
        rho_bkg_thermalandjet = 0.;
        jb_reco.reset();
        resid_rhoA = 0.;
        jb_truth.reset();
        tree->Fill();
        continue; 
    } 

    jb_truth.fill(js_reader.truth_jet, zero, true, js_reader.sigma);


    for (int n_bkg = 0; n_bkg < N_BKGSAMPLES; ++n_bkg)
    {

      auto part_M = bkgmaker(); // vector<fastjet::PseudoJet> // part_M are indiced as -1

      //background pre-Jetscape insert
      rho_bkg_thermal = bkg_est(part_M);
      auto &part_P = js_reader.all_constituents;

      //background with Jetscape insert
      part_M.insert(part_M.end(), part_P.begin(), part_P.end());
      rho_bkg_thermalandjet = bkg_est_wJetScape(part_M);

      //reconstructed ("measured" jets)
      auto jets_M = clusterer_area(part_M);

      // jet the matched jet
      select_circle.set_reference(js_reader.truth_jet);
      auto closest_jets = select_circle(jets_M);
      matched_TtoR = (closest_jets.size() != 0);
      if (!matched_TtoR) {
           jb_reco.reset();
           dR_TtoR = -1.;
           resid_rhoA = 0.;

        /* std::cout << Form("Truth jet: pt/eta/hpi %.2f/%.3f/%.2f not matched (rho is %.2f)", */ 
            /* js_reader.truth_jet.pt(), js_reader.truth_jet.eta(), js_reader.truth_jet.phi(), rho_bkg_thermalandjet) << std::endl; */
        tree->Fill();
        continue;
      }
      int i_match = 0;
      float pt_match = closest_jets[0].pt();
      for (int i = 1; i < closest_jets.size(); ++i)
      {
        if (closest_jets[i].pt() > pt_match)
        {
          i_match = i;
          pt_match = closest_jets[i].pt();
        }
      }
      auto& reco_jet = closest_jets[i_match];
      dR_TtoR = js_reader.truth_jet.delta_R(reco_jet);

      // fill in the leading jet pair
      float zero;
      jb_truth.fill(js_reader.truth_jet, zero, true, js_reader.sigma);
      jb_reco.fill(reco_jet, zero, true);

      pt_IP = js_reader.IP_jet.pt();
      eta_IP = js_reader.IP_jet.eta();
      phi_IP = js_reader.IP_jet.phi();

      resid_rhoA = (reco_jet.perp() - reco_jet.area()*rho_bkg_thermal) - js_reader.IP_jet.perp();
      tree->Fill();
    }
  }

  fout->cd();
  fout->Write();
  fout->Save();
  fout->Close();

  return 0;
};
