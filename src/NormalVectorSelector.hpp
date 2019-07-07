//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_NORMALVECTORSELECTOR_HPP
#define TINKER_NORMALVECTORSELECTOR_HPP

#include <list>
#include "VectorSelector.hpp"
#include "atom.hpp"


class NormalVectorSelector : public VectorSelector {
public:
    int initialize(const std::shared_ptr<Frame> &frame) override;

    std::vector<std::tuple<double, double, double>> calculateVectors(const std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    std::tuple<double, double, double>
    calculateVector(const std::shared_ptr<Molecule> &mol, const std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    static const std::string title() { return "Plane normal vector (define by three atoms in same molecule) selector"; }

protected:
    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;
    Atom::AtomIndenter ids3;

    std::list<std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;


    std::tuple<double, double, double> calVector(
            const std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
            const std::shared_ptr<Frame> &frame);

};


#endif //TINKER_NORMALVECTORSELECTOR_HPP
