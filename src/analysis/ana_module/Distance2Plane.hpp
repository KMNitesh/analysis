

#ifndef TINKER_DISTANCE2PLANE_HPP
#define TINKER_DISTANCE2PLANE_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "utils/PBCUtils.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class Distance2Plane : public AbstractAnalysis {
public:
    Distance2Plane();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Distance to plane"; }

private:
    std::array<AmberMask, 3> plane_marks;
    std::array<std::shared_ptr<Atom>, 3> plane_atoms;

    std::vector<AmberMask> outplane_marks;
    std::vector<std::shared_ptr<Atom>> outplane_atoms;

    PBCUtils::MolPair mol;

    struct Point {
        Point(std::tuple<double, double, double> value)
            : x(std::get<0>(value)), y(std::get<1>(value)), z(std::get<2>(value)) {}

        double x, y, z;
    };

    std::tuple<double, double, double, double> get_panel(Point p1, Point p2, Point p3);

    static double dis_pt2panel(Point pt, std::tuple<double, double, double, double> parameter);

    std::deque<double> distances;

    boost::accumulators::accumulator_set<
        double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>>
        acc;
};

#endif  // TINKER_DISTANCE2PLANE_HPP
