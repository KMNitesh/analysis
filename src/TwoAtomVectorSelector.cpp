//
// Created by xiamr on 7/5/19.
//

#include <iostream>
#include "TwoAtomVectorSelector.hpp"
#include "ThrowAssert.hpp"
#include "frame.hpp"
#include "common.hpp"

using namespace std;

int TwoAtomVectorSelector::initialize(const std::shared_ptr<Frame> &frame) {
    for (auto &mol : frame->molecule_list) {
        shared_ptr<Atom> atom1, atom2;
        for (auto &atom : mol->atom_list) {
            if (Atom::is_match(atom, ids1)) {
                atom1 = atom;
            } else if (Atom::is_match(atom, ids2)) {
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

std::vector<std::tuple<double, double, double>>
TwoAtomVectorSelector::calcaulteVectors(const std::shared_ptr<Frame> &frame) {
    std::vector<std::tuple<double, double, double>> vectors;
    for (auto &pair : pairs) {
        vectors.push_back(calVector(pair, frame));
    }
    return vectors;
}

void TwoAtomVectorSelector::readInfo() {
    std::cout << "# " << title() << " #\n";
    Atom::select1group(ids1, "Please Enter for Atom1 > ");
    Atom::select1group(ids2, "Please Enter for Atom2 > ");
}

void TwoAtomVectorSelector::print(std::ostream &os) {
    os << "# " << title() << "#\n";
    os << " Group1 > " << ids1 << "\n Group2 > " << ids2 << '\n';
}

std::tuple<double, double, double>
TwoAtomVectorSelector::calVector(std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
                                 const std::shared_ptr<Frame> &frame) {
    return get<1>(atoms)->getCoordinate() - get<0>(atoms)->getCoordinate();

}
