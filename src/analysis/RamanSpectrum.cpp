//
// Created by xiamr on 8/19/19.
//

#include <tbb/tbb.h>
#include <boost/range/algorithm.hpp>
#include "RamanSpectrum.hpp"
#include "common.hpp"
#include "atom.hpp"
#include "IRSpectrum.hpp"

RamanSpectrum::RamanSpectrum() {
    enable_tbb = true;
}

void RamanSpectrum::calculateSpectrum(const std::string &out) {

    auto time_increment_ps = choose(0.0, 100.0, "time_increment_ps [0.001 ps] :", Default(0.001));
    std::string file = choose_file("Enter Polarizable Tensor Evolution Data File : ").isExist(true);

    // Pxx, Pyy, Pzz
    std::array<std::deque<double>, 3> polarTensor;

    std::ifstream ifstream(file);

    std::array<double, 3> tensor{};
    while (ifstream) {
        ifstream >> tensor[0] >> tensor[1] >> tensor[2];
        if (polarTensor[0].empty()) {
            for (std::size_t i = 0; i < polarTensor.size(); ++i) {
                polarTensor[i].emplace_back(tensor[i]);
            }
        } else {
            for (std::size_t i = 0; i < polarTensor.size(); ++i) {
                auto last_it = --polarTensor[i].end();
                *last_it = tensor[i] - *last_it;
                polarTensor[i].emplace_back(tensor[i]);
            }
        }
    }

    for (auto &component : polarTensor) {
        component.pop_back();
    }

    std::ofstream os(out);
    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << std::string(50, '#') << '\n';
    std::array<std::vector<long double>, 3> acf;

    tbb::parallel_invoke([&acf, &polarTensor] { acf[0] = IRSpectrum::calculateAcf(polarTensor[0]); },
                         [&acf, &polarTensor] { acf[1] = IRSpectrum::calculateAcf(polarTensor[1]); },
                         [&acf, &polarTensor] { acf[2] = IRSpectrum::calculateAcf(polarTensor[2]); });


    os << std::string(50, '#') << '\n';
    os << "# " << title() << '\n';
    os << "# time_increment_ps > " << time_increment_ps << '\n';
    os << std::string(50, '#') << '\n';

    os << boost::format("%15s %15s %15s %15s\n") % "Time(ps)" % "ACF(Pxx)" % "ACF(Pyy)" % "ACF(Pzz)";

    for (std::size_t i = 0; i < acf[0].size(); ++i) {
        os << boost::format("%15.4f %15.5f %15.5f %15.5f\n")
              % (time_increment_ps * i) % acf[0][i] % acf[1][i] % acf[2][i];
    }

    os << std::string(50, '#') << '\n';


    std::array<std::vector<double>, 3> intense;

    for (std::size_t i = 0; i < polarTensor.size(); ++i) {
        intense[i] = IRSpectrum::calculateIntense(acf[i], time_increment_ps);
    }


    os << boost::format("%15s %15s %15s %15s\n")
          % "Frequency (cm-1)" % "Intensity(Pxx)" % "Intensity(Pyy)" % "Intensity(Pzz)";

    for (std::size_t i = 0; i < intense[0].size(); ++i) {
        os << boost::format("%15.3f %15.5f %15.5f %15.5f\n")
              % (i + 1) % intense[0][i] % intense[1][i] % intense[2][i];
    }

    os << std::string(50, '#') << '\n';
}