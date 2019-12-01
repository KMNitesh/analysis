//
// Created by xiamr on 8/17/19.
//

#include "IRSpectrumDeltaDipole.hpp"
#include "data_structure/frame.hpp"
#include "utils/common.hpp"

void IRSpectrumDeltaDipole::process(std::shared_ptr<Frame> &frame) {
    if (dipole_evolution.empty()) {
        dipole_evolution.emplace_back(getDipole(frame));
    } else {
        auto last_it = --dipole_evolution.end();
        auto dipole = getDipole(frame);
        *last_it = dipole - *last_it;
        dipole_evolution.emplace_back(dipole);
    }
}

void IRSpectrumDeltaDipole::print(std::ostream &os) {
    dipole_evolution.pop_back();
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << std::string(50, '#') << '\n';
    IRSpectrum::print(os);
}


void IRSpectrumDeltaDipole::calculateSpectrum(const std::string &out) {
    auto time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.001 ps] :", Default(0.001));
    std::string file = choose_file("Enter Dipole Evolution Data File : ").isExist(true);

    std::deque<std::tuple<double, double, double>> dipole_evolution;

    std::ifstream ifstream(file);

    double dplx, dply, dplz;
    while (ifstream) {
        ifstream >> dplx >> dply >> dplz;
        if (dipole_evolution.empty()) {
            dipole_evolution.emplace_back(dplx, dply, dplz);
        } else {
            auto last_it = --dipole_evolution.end();
            *last_it = std::make_tuple(dplx, dply, dplz) - *last_it;
            dipole_evolution.emplace_back(dplx, dply, dplz);
        }
    }
    dipole_evolution.pop_back();

    std::ofstream os(out);
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << std::string(50, '#') << '\n';

    printData(os, dipole_evolution, time_increment_ps);
}
