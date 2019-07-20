//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP
#define TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "BasicAnalysis.hpp"
#include "common.hpp"
#include "atom.hpp"


class Frame;

class DipoleAngleWithDistanceRange : public BasicAnalysis {
public:
    DipoleAngleWithDistanceRange() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Dipole Angle Distribution with cutoff"; }

    ~DipoleAngleWithDistanceRange() override = default;


protected:

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;


    double angle_width;

    int angle_bins;

    std::unordered_map<int, size_t> hist;

    double cutoff1, cutoff2;

    void printData(std::ostream &os) const;
};


#endif //TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP
