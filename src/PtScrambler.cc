#include "PtScrambler.h"

#include <numeric>
#include <iostream>

void PtScrambler::operator()() {
    if (opt == ScrOpt::None) {
        return;
    } else if (opt == ScrOpt::Random) {
        for (int i=0; i<N-1; ++i) {
            pos[i] = rnd.Uniform();
        }
        std::sort(pos.begin(), pos.end());
        vrat[0] = pos[0];
        for (int i=1; i<N; ++i) {
            vrat[i] = pos[i] - pos[i-1];
        }
        std::sort(vrat.rbegin(),vrat.rend());

        if (print) {
            for (int i=1; i<N; ++i) std::cout << "pt_before("<<i<<") " << rdata[i] << std::endl;
        }

        auto sum = std::accumulate(rdata.begin(), rdata.end(), 0.);
        for (int i=0; i<N; ++i) {
            rdata[i] = sum * vrat[i];
            if (print) std::cout << "pt_after("<<i<<") " << rdata[i] << std::endl;
        }
        return;
    } else if (opt == ScrOpt::Equal) {
        auto sum = std::accumulate(rdata.begin(), rdata.end(), 0.);
        auto val = sum/N;
        sum /= N;
        for (int i=0; i<N; ++i) {
            rdata[i] = val;
            if (print) std::cout << "pt_after("<<i<<") " << rdata[i] << std::endl;
        }
        return;
    } else if (opt == ScrOpt::Remove) {
        for (int i=0; i<N; ++i) {
            rdata[i] = 0.;
            if (print) std::cout << "pt_after("<<i<<") " << rdata[i] << std::endl;
        }
    }
}
