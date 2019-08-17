//
// Created by xiamr on 8/17/19.
//

#include "IRSpectrumDeltaDipole.hpp"
#include "frame.hpp"
#include "common.hpp"

void IRSpectrumDeltaDipole::process(std::shared_ptr<Frame> &frame) {
    if (dipole_evolution.empty()) {
        dipole_evolution.emplace_back(frame->getDipole());
    } else {
        auto last_it = --dipole_evolution.end();
        auto dipole = frame->getDipole();
        *last_it = dipole - *last_it;
        dipole_evolution.emplace_back(dipole);
    }
}

void IRSpectrumDeltaDipole::print(std::ostream &os) {
    dipole_evolution.pop_back();
    IRSpectrum::print(os);
}
