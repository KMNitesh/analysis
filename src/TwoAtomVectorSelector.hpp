//
// Created by xiamr on 7/5/19.
//

#ifndef TINKER_TWOATOMVECTORSELECTOR_HPP
#define TINKER_TWOATOMVECTORSELECTOR_HPP

#include <list>
#include "VectorSelector.hpp"
#include "atom.hpp"
#include "TwoAtomSelectorGrammar.hpp"


class TwoAtomVectorSelector : public VectorSelector {
public:
    int initialize(const std::shared_ptr<Frame> &frame) override;

    std::vector<std::tuple<double, double, double>> calculateVectors(const std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    void readAST(const TwoAtomVectorSelectorNode &ast);

    void print(std::ostream &os) override;

    static const std::string title() { return "Two Atom vector (define by two atoms in same molecule) selector"; }

protected:
    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;
public:
    std::tuple<double, double, double>
    calculateVector(const std::shared_ptr<Molecule> &mol, const std::shared_ptr<Frame> &frame) override;

protected:

    std::list<std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;


    std::tuple<double, double, double> calVector(
            const std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>> &atoms,
            const std::shared_ptr<Frame> &frame);
};


#endif //TINKER_TWOATOMVECTORSELECTOR_HPP
