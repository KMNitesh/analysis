//
// Created by xiamr on 6/30/19.
//

#ifndef TINKER_EQUATORIALANGLE_HPP
#define TINKER_EQUATORIALANGLE_HPP

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/Histogram.hpp"

class Frame;

class EquatorialAngle : public AbstractAnalysis {
public:
    EquatorialAngle();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Over Plane Angle Distribution with cutoff"; }

protected:
    AmberMask ids1, ids2, ids3;

    std::unordered_set<std::shared_ptr<Atom>> group1, group2, group3;

    Histogram hist;

    double cutoff1, cutoff2;

    void printData(std::ostream &os) const;
};

#endif // TINKER_EQUATORIALANGLE_HPP
