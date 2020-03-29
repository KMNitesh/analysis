
#ifndef TINKER_BONDENERGYCALCULATOR_HPP
#define TINKER_BONDENERGYCALCULATOR_HPP

#include "data_structure/atom.hpp"
#include "data_structure/frame.hpp"
#include <nlohmann/json.hpp>
#include <string_view>

class BondEnergyCalculator {
public:
    BondEnergyCalculator(const std::vector<std::shared_ptr<Atom>> &atoms, const std::shared_ptr<Frame> &frame);

    BondEnergyCalculator(const AmberMask &mask, const std::shared_ptr<Frame> &frame);

    [[nodiscard]] static std::string_view title() { return "Energy Calculator for Bonded Terms"; }

    struct Term {
        double bond, angle, dihedral, improper;
        double total() const { return bond + angle + dihedral + improper; }
    };

    Term energy(const std::shared_ptr<Frame> &frame);

    std::map<int, Term> energy_with_residue_rank(const std::shared_ptr<Frame> &frame);

    static Term energy(const AmberMask &mask, const std::shared_ptr<Frame> &frame) {
        BondEnergyCalculator calculator(mask, frame);
        return calculator.energy(frame);
    }

private:
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 2>, Frame::harmonic> *> bonds;
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 3>, Frame::harmonic> *> angles;
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 4>, Frame::pdihs> *> dihedrals;
    std::vector<const std::pair<const std::array<std::shared_ptr<Atom>, 4>, Frame::pdihs> *> improper_dihedrals;
};

inline void to_json(nlohmann::json &j, const BondEnergyCalculator::Term &term) {
    j = nlohmann::json{term.bond, term.angle, term.dihedral, term.improper};
}

inline void from_json(const nlohmann::json &j, BondEnergyCalculator::Term &term) {
    term.bond = j[0];
    term.angle = j[1];
    term.dihedral = j[2];
    term.improper = j[3];
}

#endif // TINKER_BONDENERGYCALCULATOR_HPP