//
// Created by xiamr on 8/19/19.
//

#ifndef TINKER_RAMANSPECTRUM_HPP
#define TINKER_RAMANSPECTRUM_HPP

#include "std.hpp"

class RamanSpectrum {
public:
    RamanSpectrum();

    static std::string title() { return "Raman radiation (IR) Spectrum"; }

    static void calculateSpectrum(const std::string &out);
};


#endif //TINKER_RAMANSPECTRUM_HPP
