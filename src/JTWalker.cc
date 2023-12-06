#include "JTWalker.h"
#include <iostream>

using std::cout;
using std::endl;

TreeSplit::TreeSplit(fastjet::PseudoJet& jet, const SD_criteria& SD) :
   trunk { jet }
{ 
  pt = jet.perp();
  eta = jet.eta();
  phi = (float) jet.phi_02pi();

  if (jet.has_constituents()) {
    auto parts = jet.constituents();
    n_constituents = (int) parts.size();
    int i = 0;
    /* cout << " PRE charge: " << charge << endl; */
    for (auto& P : parts) {
      charge += P.user_index();
    }
  } else {
    n_constituents = 1;
    charge = jet.user_index();
  }

  pass_pt = jet.perp() >= SD.min_pt;
  if (!pass_pt) return;

  has_parents = trunk.has_parents(branch_0, branch_1);
  if (!has_parents) return;

  if (branch_1.perp() > branch_0.perp()) {
    auto _ = branch_1;
    branch_1 = branch_0;
    branch_0 = _;
  }

  pt_0 = branch_0.perp();
  eta_0 = branch_0.eta();
  phi_0 = branch_0.phi();

  pt_1 = branch_1.perp();
  eta_1 = branch_1.eta();
  phi_1 = branch_1.phi();

  float dphi = fabs(phi_1 - phi_0);
  float deta = eta_1-eta_0;
  while (dphi > M_PI) dphi = fabs(dphi-2*M_PI);
  delta_R = pow(dphi*dphi+deta*deta,0.5);

  zg = std::min(pt_0,pt_1)/(pt_0+pt_1);
  pass_SD = (SD.beta == 0) ? (zg > SD.zg) : (zg > SD.zg * pow(delta_R/SD.R0, SD.beta));
}

TreeSplitPrinter::TreeSplitPrinter(fastjet::PseudoJet& _jet, const SD_criteria& _SD, std::string _pre_string, std::ostream& _oss) 
    : jet  {_jet}
    , SD   {_SD}
    , oss  {_oss}
{
  pre_string << _pre_string;

  TreeSplit idat { jet, SD };
  line_0 << pre_string.str() << Form("│ %-4.1f@%5.2f,%4.2f ", idat.pt, idat.eta, idat.phi);
  line_1 << pre_string.str() << Form("│ Q%3i n(%2i)      ",   (int)idat.charge, idat.n_constituents);
  line_arrow  << pre_string.str() << "└─────────────────";
  line_splits << pre_string.str() << " " << blank_string;
  line_vert  << pre_string.str() << " ";

  bool do_loop = idat.has_parents && idat.pass_pt;


  std::vector<TreeSplit> vdat {{ idat }};
  std::stack<std::string> down_streams {};

  while (vdat.back().has_parents && vdat.back().pass_pt) {
    auto& dat = vdat.back();

    TreeSplit dat1 { dat.branch_1, SD };

    std::ostringstream down_stream;
  
    bool cap_1 = !dat1.pass_pt;
    if (cap_1) {
      line_splits << Form("▼ z%4.2f ΔR(%4.2f)  ", dat.zg, dat.delta_R);
      line_vert   << Form("                 ");
    } else {
      line_splits << Form("│ z%4.2f ΔR(%4.2f)  ", dat.zg, dat.delta_R);
      line_vert   << Form("                 ");

      TreeSplitPrinter  print{dat.branch_1, SD, line_vert.str(), down_stream};
      down_streams.push(down_stream.str());
    }

    vdat.push_back( {dat.branch_0, SD} );
    auto& ndat = vdat.back();


    line_0 << Form("  %-4.1f@%5.2f,%4.2f ", ndat.pt, ndat.eta, ndat.phi);
    line_1 << Form("  Q%3i n(%2i)      ",   (int)ndat.charge, ndat.n_constituents);
  /* cout << " HOW MANY " << ndat.n_constituents << endl; */
  /* cout << " 222 MANY " << ne_-11 << Form(" Q%3i (%2i)       ",   ndat.charge, ndat.n_constituents); */
    
    line_arrow << "┬─────────────────";
    if (cap_1) line_vert << " ";
    else       line_vert << "│";

  }
  line_arrow << "─▶";
  /* print(); */
  /* oss << pre_string.str() << "│" << std::endl; */
  oss <<  line_0.str() << std::endl;
  oss <<  line_1.str() << std::endl;
  oss <<  line_arrow.str() << std::endl;
  oss <<  line_splits.str() << std::endl;
  oss <<  line_vert.str() << std::endl;
  
  while (down_streams.size() > 0) {
    oss << down_streams.top();
    down_streams.pop();
  }
  /* oss << down_streams.back().str(); */


}

/* void TreeSplitPrinter::print() { */
/*   oss << " NOT HERE " << endl; */
/*   oss << pre_string.str() << "│" << std::endl; */
/*   oss << line_0.str() << std::endl; */
/*   oss << line_1.str() << std::endl; */
/*   oss << line_arrow.str() << std::endl; */
/*   /1* if (print_branches.size() != 0) oss << line_splits.str() << std::endl; *1/ */

/*   /1* while(print_branches.size() != 0) { *1/ */
/*   /1*   print_branches.top().print(); *1/ */
/*   /1*   print_branches.pop(); *1/ */
/*   /1* } *1/ */
/* } */
