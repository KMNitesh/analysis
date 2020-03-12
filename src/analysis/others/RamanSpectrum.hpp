//
// Created by xiamr on 8/19/19.
//

#ifndef TINKER_RAMANSPECTRUM_HPP
#define TINKER_RAMANSPECTRUM_HPP

#include "utils/std.hpp"

class RamanSpectrum {
public:
    RamanSpectrum();

    [[nodiscard]] static std::string_view title() { return "Raman radiation (IR) Spectrum"; }

    static void calculateSpectrum(const std::string &out);
};

#endif  // TINKER_RAMANSPECTRUM_HPP
