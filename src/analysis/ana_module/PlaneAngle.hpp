//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_PLANEANGLE_HPP
#define TINKER_PLANEANGLE_HPP

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"

class Frame;

class PlaneAngle : public AbstractAnalysis {
public:

    PlaneAngle();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Plane Angle Distribution with cutoff"; }

protected:

    Atom::AmberMask ids1, ids2, ids3, ids4;

    std::unordered_set<std::shared_ptr<Atom>> group1, group2, group3, group4;

    std::list<std::tuple<std::shared_ptr<Atom>, std::shared_ptr<Atom>, std::shared_ptr<Atom>>> pairs;

    double angle_width;

    int angle_bins;

    std::unordered_map<int, std::size_t> hist;

    double cutoff1, cutoff2;

    void printData(std::ostream &os) const;
};


#endif //TINKER_PLANEANGLE_HPP
