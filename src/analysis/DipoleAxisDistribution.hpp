//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_DIPOLEAXISDISTRIBUTION_HPP
#define TINKER_DIPOLEAXISDISTRIBUTION_HPP

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <list>
#include "BasicAnalysis.hpp"
#include "common.hpp"
#include "atom.hpp"
#include "Histrogram.hpp"

class Frame;

class DipoleAxisDistribution : public BasicAnalysis {

public:
    DipoleAxisDistribution() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Dipole Angle Distribution with Axis"; }

protected:

    Atom::AtomIndenter ids;

    std::unordered_set<std::shared_ptr<Molecule>> group;

    Histrogram hist;

    void printData(std::ostream &os) const;

    std::tuple<double, double, double> axis_vector;

    const std::tuple<double, double, double> xaxis_vector{1, 0, 0}, yaxis_vector{0, 1, 0}, zaxis_vector{0, 0, 1};

    int axis_num;

};


#endif //TINKER_DIPOLEAXISDISTRIBUTION_HPP
