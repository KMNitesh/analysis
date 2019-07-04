//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_DIPOLEAXISDISTRIBUTION_HPP
#define TINKER_DIPOLEAXISDISTRIBUTION_HPP

#include <unordered_map>
#include <unordered_set>
#include <memory>
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

    std::unordered_set<std::shared_ptr<Atom>> group;

    Histrogram hist;

    void printData(std::ostream &os) const;

    double xr = 0.0, yr = 0.0, zr = 0.0;
};


#endif //TINKER_DIPOLEAXISDISTRIBUTION_HPP
