#ifndef PtScrambler__h
#define PtScrambler__h

#include <vector>
#include "TRandom3.h"
#include <iostream>

using std::vector;

struct PtScrambler {
    size_t N;
    vector<float>& rdata; // data to scramble
    vector<float> pos;
    vector<float> vrat;
    void operator()();
    TRandom3& rnd;
    enum ScrOpt { None = 0, Random = 1, Equal = 2, Remove = 3 };
    bool print;
    ScrOpt opt { None };


    PtScrambler(vector<float>& pt_list, TRandom3& _rand, ScrOpt _opt=ScrOpt::None, bool _print=false) 
        : N { pt_list.size() }
        , rdata { pt_list }
        , pos ( N, 0. )
        , vrat ( N, 0. )
        , rnd { _rand }
        , opt { _opt }
        , print { _print }
    { 
        pos[N-1] = 1.;
        /* std::cout << " _opt " << _opt << std::endl; */
        /* std::cout << " (IS !!) " << (_opt == ScrOpt::Random) << std::endl; */
    };
};

#endif
