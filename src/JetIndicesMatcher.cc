#include "JetIndicesMatcher.h"
#include "cmath"
#include "fastjet/PseudoJet.hh"

using namespace std;
using fastjet::PseudoJet;

JetIndicesMatcher::JetIndicesMatcher(float R) : R2 {R*R}
{}

array<unsigned int,3> JetIndicesMatcher::operator()
  ( vector<PseudoJet>& jets_T
  , vector<PseudoJet>& jets_R ) 
{

    for (unsigned int i=0; i<jets_T.size(); ++i) {
        if (i>0) {
            float diff = (jets_T[i].perp() - jets_T[i-1].perp());
        if (diff>0) {
            cout << "TRUE diff: " << diff << endl;
            cout << "-i (" << (i-1) <<")" << jets_T[i-1].perp() << " -> " << 
            " (" << i <<")" << jets_T[i].perp() << endl;
        }
        }
    }
    for (unsigned int i=0; i<jets_R.size(); ++i) {
        if (i>0) {
            float diff = (jets_R[i].perp() - jets_R[i-1].perp());
        if (diff>0) {
            cout << "MEAS diff: " << diff << endl;
            cout << "-i (" << (i-1) <<")" << jets_R[i-1].perp() << " -> " << 
            " (" << i <<")" << jets_R[i].perp() << endl;
        }
        }
    }

    i_fake.clear();
    i_miss.clear();
    i_matched.clear();

    vector<bool> is_matched (jets_R.size(),false);

    for (unsigned int T=0;T<jets_T.size();++T) {
        bool found_match { false };
        for (unsigned int R=0;R<jets_R.size();++R) {
          if (is_matched[R]) continue;
          if (jets_T[T].squared_distance(jets_R[R]) <= R2) {
            found_match = true;
            is_matched[R] = true;
            i_matched.push_back({T,R});
            break;
          }
        }
        if (!found_match) i_miss.push_back(T);
    }
    for (unsigned int R=0;R<jets_R.size();++R) {
        if (!is_matched[R]) i_fake.push_back(R);
    }
    return {static_cast<unsigned int>(i_matched.size()),
            static_cast<unsigned int>(i_miss.size()),
            static_cast<unsigned int>(i_fake.size())
    };
}
