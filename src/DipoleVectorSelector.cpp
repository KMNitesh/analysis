//
// Created by xiamr on 7/5/19.
//
#include "frame.hpp"
#include "DipoleVectorSelector.hpp"
#include "ThrowAssert.hpp"

using namespace std;

int DipoleVectorSelector::initialize(const std::shared_ptr<Frame> &frame) {
    if (!selected_mols.empty()) return selected_mols.size();
    for (auto &atom : frame->atom_list) {
        if (Atom::is_match(atom, id)) {
            selected_mols.insert(atom->molecule.lock());
        }
    }
    throw_assert(!selected_mols.empty(), "Atom selection semtatic error");
    return selected_mols.size();
}

std::vector<std::tuple<double, double, double>>
DipoleVectorSelector::calculateVectors(const std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> vectors;

    for (auto &mol : selected_mols) {
        vectors.push_back(calculateVector(mol, frame));
    }
    return vectors;
}

void DipoleVectorSelector::readInfo() {
    std::cout << "# " << title() << " #\n";
    Atom::select1group(id, "Please Enter for Atom > ");
}

void DipoleVectorSelector::print(std::ostream &os) {
    os << "# " << title() << "#\n";
    os << " Group1 > " << id << '\n';
}

tuple<double, double, double>
DipoleVectorSelector::calculateVector(const std::shared_ptr<Molecule> &mol, const std::shared_ptr<Frame> &frame) {
    auto r = mol->calc_dipole(frame);
    r /= vector_norm(r);
    return r;
}

void DipoleVectorSelector::setParameters(const Atom::Node &id) {
    this->id = id;
}
