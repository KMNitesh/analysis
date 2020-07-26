
#ifndef TINKER_BONDENERGYCOMPARE_HPP
#define TINKER_BONDENERGYCOMPARE_HPP

#include "utils/std.hpp"

class BondEnergyCompare {
public:
    [[nodiscard]] static std::string_view title() { return "Bond Energy Compare"; }

    static void process();
};


#endif // TINKER_BONDENERGYCOMPARE_HPP