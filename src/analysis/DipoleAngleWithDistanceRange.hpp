//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP
#define TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "AbstractAnalysis.hpp"
#include "atom.hpp"


class Frame;

class DipoleAngleWithDistanceRange : public AbstractAnalysis {
public:

    DipoleAngleWithDistanceRange();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Dipole Angle Distribution with cutoff"; }

protected:

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;


    double angle_width;

    int angle_bins;

    std::unordered_map<int, size_t> hist;

    double cutoff1, cutoff2;

    void printData(std::ostream &os) const;
};


#endif //TINKER_DIPOLEANGLEWITHDISTANCERANGE_HPP
