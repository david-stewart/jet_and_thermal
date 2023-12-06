#ifndef JTWalker_h
#define JTWalker_h

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"


#include "TString.h"
#include <iostream>
#include <sstream>
#include <stack>

struct SD_criteria {
  float zg      ;
  float beta    ;
  float R0      ;
  float min_pt  ;
  SD_criteria ( float _zg=0.3, float _beta=0., float _R0=0.4, float _min_pt=1. ) :
    zg{_zg}
  , beta{_beta}
  , R0 {_R0}
  , min_pt { _min_pt}
  {}
};

struct TreeSplit {
  fastjet::PseudoJet& trunk; // 0
  float pt, eta, phi;
  bool pass_pt;

  bool has_parents = false;
  bool pass_SD = false; 

  float pt_0 {}, eta_0 {}, phi_0 {};
  float pt_1 {}, eta_1 {}, phi_1 {};
  float delta_R;
  float zg;

  int n_constituents;
  float charge {0};

  fastjet::PseudoJet branch_0{};
  fastjet::PseudoJet branch_1{};

  TreeSplit(fastjet::PseudoJet& jet, const SD_criteria& criteria);
};

struct TreeSplitPrinter {
  fastjet::PseudoJet& jet;
  const SD_criteria& SD;
  std::ostream& oss;

  std::ostringstream pre_string {};

  std::ostringstream line_0 {}; // pt eta phi line
  std::ostringstream line_1 {}; // sumch Nconstituents line
  std::ostringstream line_arrow  {};  // 
  std::ostringstream line_vert   {};  // 
  std::ostringstream line_splits {}; // branching stats

  std::stack<TreeSplitPrinter> print_branches {};

  /* inline static const std::string hh_string    = "─────────────────"; */
  inline static const std::string blank_string = "                 ";

  /* std::string ptetaphi_string();    // return pT@eta,phi { return Form(" %4.1f@%4.2F,%4.2f ", walker.current_jet.perp(), walker.current_jet. */
  /* std::string sumch_nperp_string(); // return pT@eta,phi { return Form(" %4.1f@%4.2F,%4.2f ", walker.current_jet.perp(), walker.current_jet. */


  TreeSplitPrinter(fastjet::PseudoJet& jet, const SD_criteria& criteria, std::string _pre_string="", std::ostream& oss=std::cout);
  void print();
};

#endif

/*

  Key: 

  pT@eta,phi  =  ##.#@#.##,#.## 
  Σch (Nch)   =  Σ## (##)

 │
 │ ##.#@##.##,#.##  ##.#@##.##,#.##  ##.#@##.##,#.##   ##.#@#.##,#.##    ##.#@#.##,#.##
 │ Σ### (##)        Σ## (##)         Σ## (##)         Σ## (##)          Σ## (##)       
 └────────────────┬────────────────┬────────────────┬────────────────┬─────────────────▶ 
                  ▼ z#.## ΔR(#.##) │ z#.## ΔR(#.##) │ z#.## ΔR(#.##) │  z#.## ΔR#.##  
                                   │                │                │               
                                   │                │                │  ##.#@#.##,#.##            
                                   │                │                │  Σ## (##)                                          
                  │                │                │                └─────────────────▶
                  │                │                │                    
                  │                │                │ ##.#@#.##,#.##  ##.#@#.##,#.## 
                  │                │                │ Σ## (##)         Σ## (##)       
                  │                │                └────────────────┬─────────────────▶ 
                  │                │                                 │  z#.## ΔR#.##     
                  │                │                                 │               
                  │                │                                 │  ##.#@#.##,#.##  ##.#@#.##,#.##      
                  │                │                                 │  Σ## (##)        Σ## (##)         
                  │                │                                 └─────────────────┬─────────────────▶
                  │                │                                                   │  z#.## ΔR#.##    
                  │                │                                                   │
                  │                │                                                   │  ##.#@#.##,#.## 
                  │                │                                                   │  Σ## (##)       
                  │                │                                                   └─────────────────▶ 
                  │                │                                  
                  │                │ ##.#@##.##,#.##   ##.#@##.##,#.##                                     
                  │                │  Σ## (##)         Σ## (##)                                           
                  │                └──────────────────┬─────────────────▶                                 
                  │                                   │  z#.## ΔR(#.##)  
                  │                                   │                 
                  │                                   │  ##.#@#.##,#.## 
                  │                                   │  Σ## (##)       
                  │                                   └─────────────────▶
                  │ 
                  │ ##.#@##.##,#.##  
                  │  Σ## (##)         Σ## (##)          
                  └──────────────────┬─────────────────▶
                                     │  z#.## ΔR(#.##)  
                                     │                 
                                     │  ##.#@#.##,#.## 
                                     │  Σ## (##)       
                                     └─────────────────▶
*/
