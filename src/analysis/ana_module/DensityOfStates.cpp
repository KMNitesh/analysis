//
// Created by xiamr on 8/13/19.
//

#include "DensityOfStates.hpp"

#include "utils/common.hpp"

void DensityOfStates::print(std::ostream &os) {
    VelocityAutocorrelationFunction::print(os);

    constexpr double lightSpeed = 2.99792458E-2;

    const double factor = 2.0 * pi * lightSpeed;

    std::vector<double> intense(300, 0);
    for (size_t f = 0; f < intense.size(); ++f) {
        auto freq = f + 1;
        for (size_t k = 0; k < acf.size(); ++k) {
            intense[f] += acf[k] * cos(freq * factor * time_increment_ps * k);
        }
    }

    auto max_intense = *std::max_element(intense.begin(), intense.end());

    std::for_each(intense.begin(), intense.end(), [max_intense](auto &i) { i /= max_intense; });

    os << boost::format("%15s %15s\n") % "Frequency (cm-1)" % "Intensity";

    for (std::size_t i = 0; i < intense.size(); ++i) {
        os << boost::format("%15.3f %15.5f\n") % (i + 1) % intense[i];
    }

    os << std::string(50, '#') << '\n';
}
