//
// Created by xiamr on 7/5/19.
//

#ifndef TINKER_DIPOLEVECTORSELECTOR_HPP
#define TINKER_DIPOLEVECTORSELECTOR_HPP

#include <set>
#include "VectorSelector.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "common.hpp"

class DipoleVectorSelector : public VectorSelector {
public:
    DipoleVectorSelector() { enable_forcefield = true; }

    int initialize(const std::shared_ptr<Frame> &frame) override;

    std::vector<std::tuple<double, double, double>> calculateVectors(const std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    void setParameters(const Atom::Node &id);

    std::tuple<double, double, double>
    calculateVector(const std::shared_ptr<Molecule> &mol, const std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    static const std::string title() { return "Dipole vector (define by molecule that has selected atom) selector"; }

protected:
    Atom::AtomIndenter id;

    std::set<std::shared_ptr<Molecule>> selected_mols;
};


#endif //TINKER_DIPOLEVECTORSELECTOR_HPP
