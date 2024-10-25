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

void add_dot_root(string& name) {
  if ( (name.size() < 5) ||
     (name.compare(name.length()-5, 5, ".root") != 0)
  ) {
    name+=".root";
  } 
  return;
}

int main(int nargs, char **argv)
{
  /*
   * arguments:
   *   1: input file
   *   2: max_events (defualts to -1)
   *   3: output file
   *   4: input_bkg file
   *   5: bool -- fill the vector of cpt
   *   6: N_bkgsamples per jet
   *   7: bkg_seed (defaults to -100 and is set to 0 -- therefore random)
   */

    /* In order to run background events, do something like:
    ./bin/IP_match_trees INPUT_PROBE 10 probe_file_out.root dat/hydro_31K.root false 1 0
    */

  // running parmaeters
  string       f_in         {(nargs > 1) ? argv[1] : "test_in.root"};
  // use "INPUT_PROBE" to get this code to run 30 GeV probes into the background
  const int    MAX_NEVENTS  {(nargs > 2) ? atoi(argv[2]) : -1};
  string       f_out        {(nargs > 3) ? argv[3] : "test_out.root"};
  string       bulk_file    {(nargs > 4) ? argv[4] : "none"}; // is_hydro
  bool         fill_vec_cpt {(nargs > 5) ? static_cast<bool>(atoi(argv[5])) : false};
  const int    N_BKGSAMPLES {(nargs > 6) ? atoi(argv[6]) : 1};
  int          bkg_seed     {(nargs > 7) ? atoi(argv[6]) : -100};  // defaults to equal pythia seed

  std::cout << " Input: f_in(" << f_in <<")" << std::endl
            <<" MAX_NEVENTS("<<MAX_NEVENTS <<")" << std::endl
            <<" f_out("<<f_out <<")" << std::endl
            <<" bulk_file("<<bulk_file <<")" << std::endl
            <<" fill_vec_cpt("<<fill_vec_cpt <<")" << std::endl
            <<" N_BKGSAMPLES("<<N_BKGSAMPLES <<")" << std::endl
            <<" bkg_seed("<<bkg_seed << ")" << std::endl; 

  bool use_probe = (f_in == "INPUT_PROBE");
  if (use_probe) { f_in = "test_in"; };

  add_dot_root(f_in);
  add_dot_root(f_out);

  TFile* bulk_tfile = nullptr;
  TTree* bulk_ttree = nullptr;
  bool   separate_bulk = false;
  if (bulk_file != "is_hydro") 
  {
    add_dot_root(bulk_file);
    std::cout << " Reading in separate hydro background file: " << bulk_file << std::endl;
    bulk_tfile = new TFile(bulk_file.c_str(), "read");
    bulk_ttree = (TTree*)bulk_tfile->Get("T");
    separate_bulk = true;
    assert(bulk_ttree != nullptr);
  } 

  cout << " Input: " << f_in << "->" << f_out <<" vec("<<fill_vec_cpt<<") nemb="<<N_BKGSAMPLES<<" ev="<<MAX_NEVENTS<<endl;

  if (bkg_seed == -100) bkg_seed = 0.;

  TFile *fout = new TFile((f_out).c_str(), "recreate");
  TTree *tree = new TTree("T", "Trees for the initiating parton");


  TH1D* dR_rap = new TH1D("dR_rap","dR using rapidity;#sqrt{#Delta#y+#Delta#phi};Truth to Reco",100, 0., 0.4);
  TH1D* dR_eta = new TH1D("dR_eta","dR using pseudo-rapidity;#sqrt{#Delta#eta+#Delta#phi};Truth to Reco",100, 0., 0.4);

  TH1D* IP_rap = new TH1D("IP_rap","IP_rap;#mathit{y};",100, -1.5, 1.5);
  TH1D* IP_eta = new TH1D("IP_eta","IP_eta;#mathit{y};",100, -1.5, 1.5);
  TH1D* IP_deltaetarap = new TH1D("IP_deltaetarap","IP eta - rap;#eta - #mathit{y};",100, -0.15, 0.15);

  TH1D* truth_rap = new TH1D("truth_rap",";truth jet rapidity;",100, -1.5, 1.5);
  TH1D* truth_eta = new TH1D("truth_eta",";truth jet pseudorapidity;",100, -1.5, 1.5);
  TH1D* truth_deltaetarap = new TH1D("truth_deltaetarap","truth eta - rap;#eta - #mathit{y};",100, -0.15, 0.15);

  TH1D* reco_rap = new TH1D("reco_rap",";reco jet rapidity;",100, -1.5, 1.5);
  TH1D* reco_eta = new TH1D("reco_eta",";reco jet pseudorapidity;",100, -1.5, 1.5);
  TH1D* reco_deltaetarap = new TH1D("reco_deltaetarap","reco eta - rap;#eta - #mathit{y};",100, -0.45, 0.45);
  TH2D* reco_etarap = new TH2D("reco_etarap","reco eta - rap;#eta;#textit{y}",100, -1.5, 1.5, 100, -1.5, 1.5);

  IP_geom_matcher js_reader(f_in, use_probe, bkg_seed, MAX_NEVENTS);
  if (bulk_ttree == nullptr) { 
    bulk_ttree = js_reader.tree; 
    separate_bulk = false; 
  }

  float pt_IP;
  float eta_IP;
  float phi_IP;

  bool matched_TtoR  = false;  // truth to reco
  bool matched_IPtoT = false; // IP to truth

  float dR_IPtoT = -1.;
  float dR_TtoR  = -1.;

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
  bkgmaker.init(bulk_ttree, separate_bulk);

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
    if (MAX_NEVENTS != -1 && js_reader.i_event > MAX_NEVENTS)
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

    IP_rap->Fill(js_reader.IP_jet.rap());
    IP_eta->Fill(js_reader.IP_jet.eta());
    IP_deltaetarap->Fill(js_reader.IP_jet.eta() -js_reader.IP_jet.rap());

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

    truth_rap->Fill(js_reader.truth_jet.rap());
    truth_eta->Fill(js_reader.truth_jet.eta());
    truth_deltaetarap->Fill(js_reader.truth_jet.eta()-js_reader.truth_jet.rap());
    /* std::cout << " " << js_reader.truth_jet.rap() << " vs " << js_reader.truth_jet.eta() << std::endl; */

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

      reco_rap->Fill(reco_jet.rap());
      reco_eta->Fill(reco_jet.eta());
      reco_deltaetarap->Fill(reco_jet.eta()-reco_jet.rap());
      reco_etarap->Fill(     reco_jet.eta(),reco_jet.rap());

      dR_TtoR = js_reader.truth_jet.delta_R(reco_jet);

      dR_rap->Fill(dR_TtoR);

      float temp_deta = js_reader.truth_jet.eta() - reco_jet.eta();
      float temp_dphi = js_reader.truth_jet.phi() - reco_jet.phi();
      while (temp_dphi > M_PI) temp_dphi -= 2*M_PI;
      while (temp_dphi < -M_PI) temp_dphi += 2*M_PI;
      float temp_dR = sqrt(temp_deta*temp_deta+temp_dphi*temp_dphi);
      dR_eta->Fill(temp_dR);
      /* std::cout << " FIXME dR_TtoR " << dR_TtoR << " " << fixme_dR << " dR_with_rap: " << fixme_dR_rap << std::endl; */



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
