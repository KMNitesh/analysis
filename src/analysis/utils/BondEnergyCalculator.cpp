
#include "BondEnergyCalculator.hpp"
#include "data_structure/frame.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"
#include <boost/range/algorithm.hpp>

BondEnergyCalculator::BondEnergyCalculator(const std::vector<std::shared_ptr<Atom>> &atoms,
                                           const std::shared_ptr<Frame> &frame) {

    std::set atoms_set(std::begin(atoms), std::end(atoms));

    for (const auto &bond : frame->f_bond_params) {
        if (atoms_set.contains(bond.first[0]) and atoms_set.contains(bond.first[1])) {
            bonds.push_back(&bond);
        }
    }

    for (const auto &angle : frame->f_angle_params) {
        if (atoms_set.contains(angle.first[0]) and atoms_set.contains(angle.first[1]) and
            atoms_set.contains(angle.first[2])) {
            angles.push_back(&angle);
        }
    }
    for (const auto &dihedral : frame->f_dihedral_params) {
        if (atoms_set.contains(dihedral.first[0]) and atoms_set.contains(dihedral.first[1]) and
            atoms_set.contains(dihedral.first[2]) and atoms_set.contains(dihedral.first[3])) {
            dihedrals.push_back(&dihedral);
        }
    }
}

BondEnergyCalculator::BondEnergyCalculator(const AmberMask &mask, const std::shared_ptr<Frame> &frame)
    : BondEnergyCalculator(PBCUtils::find_atoms(mask, frame), frame) {}

BondEnergyCalculator::Term BondEnergyCalculator::energy(const std::shared_ptr<Frame> &frame) {

    Term term{};
    for (auto it = bonds.begin(); it != bonds.end(); ++it) {
        const auto [atoms, param] = **it;
        auto r = atom_distance(atoms, frame) - param.rA;
        term.bond += 0.5 * param.krA * r * r;
    }

    for (auto it = angles.begin(); it != angles.end(); ++it) {
        const auto [atoms, param] = **it;
        auto theta = (atom_angle(atoms, frame) - param.rA) / radian;
        term.angle += 0.5 * param.krA * theta * theta;
    }
    for (auto it = dihedrals.begin(); it != dihedrals.end(); ++it) {
        const auto [atoms, param] = **it;
        auto cos = std::cos((atom_dihedral(atoms, frame) * param.mult - param.phiA) / radian);
        term.dihedral += param.cpA * (1 + cos);
    }
    return term;
}