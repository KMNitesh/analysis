
#include "BondEnergyCalculator.hpp"
#include "data_structure/frame.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"
#include <boost/range/algorithm.hpp>

BondEnergyCalculator::Term BondEnergyCalculator::energy(const AmberMask &mask, const std::shared_ptr<Frame> &frame) {
    return energy(PBCUtils::find_atoms(mask, frame), frame);
}

BondEnergyCalculator::Term BondEnergyCalculator::energy(const std::vector<std::shared_ptr<Atom>> &atoms,
                                                        const std::shared_ptr<Frame> &frame) {
    std::set atoms_set(std::begin(atoms), std::end(atoms));

    Term term{};
    for (const auto &[key, param] : frame->f_bond_params) {
        std::array atom{frame->atom_map[key[0]], frame->atom_map[key[1]]};
        if (atoms_set.contains(atom[0]) and atoms_set.contains(atom[1])) {
            auto r = atom_distance(atom, frame) - param.rA;
            term.bond += 0.5 * param.krA * r * r;
        }
    }

    for (const auto &[key, param] : frame->f_angle_params) {
        std::array atom{frame->atom_map[key[0]], frame->atom_map[key[1]], frame->atom_map[key[2]]};
        if (atoms_set.contains(atom[0]) and atoms_set.contains(atom[1]) and atoms_set.contains(atom[2])) {
            auto theta = (atom_angle(atom, frame) - param.rA) / radian;
            term.angle += 0.5 * param.krA * theta * theta;
        }
    }
    for (const auto &[key, param] : frame->f_dihedral_params) {
        std::array atom{frame->atom_map[key[0]], frame->atom_map[key[1]], frame->atom_map[key[2]],
                        frame->atom_map[key[3]]};
        if (atoms_set.contains(atom[0]) and atoms_set.contains(atom[1]) and atoms_set.contains(atom[2]) and
            atoms_set.contains(atom[3])) {
            auto cos = std::cos((atom_dihedral(atom, frame) * param.mult - param.phiA) / radian);
            term.dihedral += param.cpA * (1 + cos);
        }
    }

    return term;
}