#ifndef JetInidicesMatcher__h
#define JetInidicesMatcher__h

#include <vector>
#include <array>
#include "fastjet/PseudoJet.hh"

struct JetIndicesMatcher {
    // a simple class that matches jets ranked on pT in eta-phi space
    // This is done by matching the first jets (already sorted by pT)
    // with phi-eta distance less than R
    const double R2; // distance in phi-eta squared

    JetIndicesMatcher(float R=0.4);

    // values calculated
    std::vector<unsigned int> i_fake{};
    std::vector<unsigned int> i_miss{};
    std::vector<std::pair < unsigned int, unsigned int> > i_matched{};

    std::array<unsigned int, 3> operator()(
    std::vector<fastjet::PseudoJet>& jets_truth,
    std::vector<fastjet::PseudoJet>& jets_wbkg);

    std::array<unsigned int, 3> operator()(
       fastjet::PseudoJet& lead_jet,
       std::vector<fastjet::PseudoJet>& jets_wbkg) 
    {
        std::vector<fastjet::PseudoJet> jets_truth { lead_jet };
        return this->operator()(jets_truth, jets_wbkg);
    };
};
#endif
