//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_NORMALVECTORSELECTOR_HPP
#define TINKER_NORMALVECTORSELECTOR_HPP

#include <list>

#include "data_structure/atom.hpp"
#include "utils/VectorSelector.hpp"

class NormalVectorSelector : public VectorSelector {
public:
    int initialize(const std::shared_ptr<Frame> &frame) override;

    std::vector<std::tuple<double, double, double>> calculateVectors(const std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    void setParameters(const AmberMask &id1, const AmberMask &id2, const AmberMask &id3);

    std::tuple<double, double, double> calculateVector(const std::shared_ptr<Molecule> &mol,
                                                       const std::shared_ptr<Frame> &frame) override;

    std::string description() override;

    void print(std::ostream &os) override;

    static const std::string title() { return "Plane normal vector (define by three atoms in same molecule) selector"; }

protected:
    AmberMask ids1;
    AmberMask ids2;
    AmberMask ids3;

    std::list<std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;

    std::tuple<double, double, double> calVector(
        const std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
        const std::shared_ptr<Frame> &frame);
};

#endif  // TINKER_NORMALVECTORSELECTOR_HPP
