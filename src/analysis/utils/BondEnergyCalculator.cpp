
#include "BondEnergyCalculator.hpp"
#include "data_structure/frame.hpp"
#include "utils/PBCUtils.hpp"
#include "utils/common.hpp"
#include <boost/range/algorithm.hpp>

namespace {

template <typename T, typename U> bool all_in_set(const T &sequecne, const std::set<U> &set) {
    for (auto &i : sequecne) {
        if (!set.contains(i))
            return false;
    }
    return true;
}
} // namespace

BondEnergyCalculator::BondEnergyCalculator(const std::vector<std::shared_ptr<Atom>> &atoms,
                                           const std::shared_ptr<Frame> &frame) {

    std::set atoms_set(std::begin(atoms), std::end(atoms));

    for (const auto &bond : frame->f_bond_params) {
        if (all_in_set(bond.first, atoms_set)) {
            bonds.push_back(std::cref(bond));
        }
    }
    for (const auto &angle : frame->f_angle_params) {
        if (all_in_set(angle.first, atoms_set)) {
            angles.push_back(std::cref(angle));
        }
    }
    for (const auto &dihedral : frame->f_dihedral_params) {
        if (all_in_set(dihedral.first, atoms_set)) {
            dihedrals.push_back(std::cref(dihedral));
        }
    }
    for (const auto &improper : frame->f_improper_dihedral_params) {
        if (all_in_set(improper.first, atoms_set)) {
            improper_dihedrals.push_back(std::cref(improper));
        }
    }
}

BondEnergyCalculator::BondEnergyCalculator(const AmberMask &mask, const std::shared_ptr<Frame> &frame)
    : BondEnergyCalculator(PBCUtils::find_atoms(mask, frame), frame) {}

BondEnergyCalculator::Term BondEnergyCalculator::energy(const std::shared_ptr<Frame> &frame) {

    Term term{};
    for (auto it = bonds.begin(); it != bonds.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto r = atom_distance(atoms, frame) - param.rA;
        term.bond += 0.5 * param.krA * r * r;
    }
    for (auto it = angles.begin(); it != angles.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto theta = (atom_angle(atoms, frame) - param.rA) * degree;
        term.angle += 0.5 * param.krA * theta * theta;
    }
    for (auto it = dihedrals.begin(); it != dihedrals.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto cos = std::cos((atom_dihedral(atoms, frame) * param.mult - param.phiA) * degree);
        term.dihedral += param.cpA * (1 + cos);
    }
    for (auto it = improper_dihedrals.begin(); it != improper_dihedrals.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto cos = std::cos((atom_dihedral(atoms, frame) * param.mult - param.phiA) * degree);
        term.improper += param.cpA * (1 + cos);
    }
    return term;
}

std::map<int, BondEnergyCalculator::Term>
BondEnergyCalculator::energy_with_residue_rank(const std::shared_ptr<Frame> &frame) {

    std::map<int, Term> terms_map;
    for (auto it = bonds.begin(); it != bonds.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto r = atom_distance(atoms, frame) - param.rA;
        auto half_bond = 0.25 * param.krA * r * r;
        for (auto &atom : atoms)
            terms_map[atom->residue_num.get()].bond += half_bond;
    }
    for (auto it = angles.begin(); it != angles.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto theta = (atom_angle(atoms, frame) - param.rA) * degree;
        auto angle = 1 / 3.0 * 0.5 * param.krA * theta * theta;
        for (auto &atom : atoms)
            terms_map[atom->residue_num.get()].angle += angle;
    }
    for (auto it = dihedrals.begin(); it != dihedrals.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto cos = std::cos((atom_dihedral(atoms, frame) * param.mult - param.phiA) * degree);
        auto dihedral = 0.25 * param.cpA * (1 + cos);
        for (auto &atom : atoms)
            terms_map[atom->residue_num.get()].dihedral += dihedral;
    }
    for (auto it = improper_dihedrals.begin(); it != improper_dihedrals.end(); ++it) {
        const auto &[atoms, param] = it->get();
        auto cos = std::cos((atom_dihedral(atoms, frame) * param.mult - param.phiA) * degree);
        auto improper = 0.25 * param.cpA * (1 + cos);
        for (auto &atom : atoms)
            terms_map[atom->residue_num.get()].improper += improper;
    }
    return terms_map;
}