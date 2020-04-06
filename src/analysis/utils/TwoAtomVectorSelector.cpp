//
// Created by xiamr on 7/5/19.
//

#include "TwoAtomVectorSelector.hpp"

#include <iostream>

#include "ThrowAssert.hpp"
#include "common.hpp"
#include "data_structure/frame.hpp"
#include "data_structure/molecule.hpp"

using namespace std;

int TwoAtomVectorSelector::initialize(const std::shared_ptr<Frame> &frame) {
    if (!pairs.empty()) return pairs.size();
    for (auto &mol : frame->molecule_list) {
        shared_ptr<Atom> atom1, atom2;
        for (auto &atom : mol->atom_list) {
            if (is_match(atom, mask1)) {
                atom1 = atom;
            } else if (is_match(atom, mask2)) {
                atom2 = atom;
            }
        }

        throw_assert((atom1 && atom2) or (!atom1 && !atom2), "Atom selection semantic error");
        if (atom1 && atom2) {
            pairs.emplace_back(atom1, atom2);
        }
    }

    throw_assert(!pairs.empty(), "Can not empty");
    return pairs.size();
}

std::vector<std::tuple<double, double, double>> TwoAtomVectorSelector::calculateVectors(
    const std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> vectors;
    for (auto &pair : pairs) {
        vectors.push_back(calVector(pair, frame));
    }
    return vectors;
}

void TwoAtomVectorSelector::readInfo() {
    std::cout << "# " << title() << " #\n";
    select1group(mask1, "Please Enter for Atom1(start) > ");
    select1group(mask2, "Please Enter for Atom2(end)   > ");
}

void TwoAtomVectorSelector::print(std::ostream &os) {
    os << "# " << title() << "#\n";
    os << "# Group1 > " << mask1 << "\n Group2 > " << mask2 << '\n';
}

std::tuple<double, double, double> TwoAtomVectorSelector::calVector(
    const std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms, const std::shared_ptr<Frame> &frame) {
    auto r = get<1>(atoms)->getCoordinate() - get<0>(atoms)->getCoordinate();
    frame->image(r);
    r /= vector_norm(r);
    return r;
}

tuple<double, double, double> TwoAtomVectorSelector::calculateVector(const std::shared_ptr<Molecule> &mol,
                                                                     const std::shared_ptr<Frame> &frame) {
    shared_ptr<Atom> atom1, atom2;
    for (auto &atom : mol->atom_list) {
        if (is_match(atom, mask1)) {
            atom1 = atom;
        } else if (is_match(atom, mask2)) {
            atom2 = atom;
        }
    }

    throw_assert(atom1 && atom2, "Atom selection semantic error");
    return calVector({atom1, atom2}, frame);
}

void TwoAtomVectorSelector::setParameters(const AmberMask &id1, const AmberMask &id2) {
    this->mask1 = id1;
    this->mask2 = id2;
}

string TwoAtomVectorSelector::description() {
    stringstream ss;
    ss << "TwoAtomVector ( [" << mask1 << "] , [" << mask2 << "] )";
    return ss.str();
}
