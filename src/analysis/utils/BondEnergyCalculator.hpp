
#ifndef TINKER_BONDENERGYCALCULATOR_HPP
#define TINKER_BONDENERGYCALCULATOR_HPP

#include <string_view>
#include "data_structure/atom.hpp"

class BondEnergyCalculator {
public:

    [[nodiscard]] static std::string_view title() { return "Energy Calculator for Bonded Terms"; }

    struct Term
    {
        double bond, angle, dihedral;
    };
    
    static Term energy(const std::vector<std::shared_ptr<Atom>> &atoms, const std::shared_ptr<Frame> &frame);

    static Term energy(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

};

#endif // TINKER_BONDENERGYCALCULATOR_HPP