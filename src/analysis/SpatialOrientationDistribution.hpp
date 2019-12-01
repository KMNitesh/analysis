//
// Created by xiamr on 8/5/19.
//

#ifndef TINKER_SPATIALORIENTATIONDISTRIBUTION_HPP
#define TINKER_SPATIALORIENTATIONDISTRIBUTION_HPP

#include <unordered_set>
#include <deque>
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;
class Molecule;

class SpatialOrientationDistribution : public AbstractAnalysis {
public:

    SpatialOrientationDistribution();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Spatial Orientation Distribution Analysis"; }

    // throw AssertionFailureException
    [[nodiscard]] static double calculatePhiAngle(const std::tuple<double, double, double> &vector);

    // throw AssertionFailureException
    [[nodiscard]] static double calculateThetaAngle(const std::tuple<double, double, double> &vector);

protected:

    Atom::AmberMask ids;

    std::unordered_set<std::shared_ptr<Molecule>> group;

    std::deque<std::tuple<double, double, double>> normal_vectors;
};


#endif //TINKER_SPATIALORIENTATIONDISTRIBUTION_HPP
