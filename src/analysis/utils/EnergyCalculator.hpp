#ifndef TINKER_ENERGYCALCULATOR_HPP
#define TINKER_ENERGYCALCULATOR_HPP

#include "utils/std.hpp"

#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class EnergyCalculator {
public:
    struct EnergyTerm {
        double ele, lj;
    };

    EnergyTerm calculate_energy(const std::shared_ptr<Frame> &frame);

    void setMask(AmberMask &mask1, AmberMask &mask2, const std::shared_ptr<Frame> &frame);

private:
    std::vector<std::shared_ptr<Atom>> group1, group2;
};

#endif // TINKER_ENERGYCALCULATOR_HPP
