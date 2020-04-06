//
// Created by xiamr on 7/5/19.
//
#include "DipoleVectorSelector.hpp"

#include "data_structure/frame.hpp"
#include "utils/ThrowAssert.hpp"
#include "utils/common.hpp"

DipoleVectorSelector::DipoleVectorSelector() { enable_forcefield = true; }

int DipoleVectorSelector::initialize(const std::shared_ptr<Frame> &frame) {
    if (!selected_mols.empty()) return selected_mols.size();
    for (auto &atom : frame->atom_list) {
        if (is_match(atom, amberMask)) {
            selected_mols.insert(atom->molecule.lock());
        }
    }
    throw_assert(!selected_mols.empty(), "Atom selection semtatic error");
    return selected_mols.size();
}

std::vector<std::tuple<double, double, double>> DipoleVectorSelector::calculateVectors(
    const std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> vectors;

    for (auto &mol : selected_mols) {
        vectors.push_back(calculateVector(mol, frame));
    }
    return vectors;
}

void DipoleVectorSelector::readInfo() {
    std::cout << "# " << title() << " #\n";
    select1group(amberMask, "Please Enter for Atom > ");
}

void DipoleVectorSelector::print(std::ostream &os) {
    os << "# " << title() << "#\n";
    os << "# Group1 > " << amberMask << '\n';
}

std::tuple<double, double, double> DipoleVectorSelector::calculateVector(const std::shared_ptr<Molecule> &mol,
                                                                         const std::shared_ptr<Frame> &frame) {
    auto r = mol->calc_dipole(frame);
    r /= vector_norm(r);
    return r;
}

void DipoleVectorSelector::setParameters(const AmberMask &id) { this->amberMask = id; }

std::string DipoleVectorSelector::description() { return "DipoleVector ( [" + to_string(amberMask) + "] )"; }
