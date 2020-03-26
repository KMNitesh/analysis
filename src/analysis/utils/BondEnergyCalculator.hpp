
#ifndef TINKER_BONDENERGYCALCULATOR_HPP
#define TINKER_BONDENERGYCALCULATOR_HPP

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include <string_view>

class BondEnergyCalculator {
public:
    BondEnergyCalculator(const std::vector<std::shared_ptr<Atom>> &atoms, const std::shared_ptr<Frame> &frame);

    BondEnergyCalculator(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

    [[nodiscard]] static std::string_view title() { return "Energy Calculator for Bonded Terms"; }

    struct Term {
        double bond, angle, dihedral;
    };

    Term energy(const std::shared_ptr<Frame> &frame);

    static Term energy(const AmberMask &mask, const std::shared_ptr<Frame> &frame) {
        BondEnergyCalculator calculator(mask, frame);
        return calculator.energy(frame);
    }

private:
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 2>, Frame::harmonic> *> bonds;
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 3>, Frame::harmonic> *> angles;
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 4>, Frame::pdihs> *> dihedrals;
};

#endif // TINKER_BONDENERGYCALCULATOR_HPP