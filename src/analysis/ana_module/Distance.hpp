//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DISTANCE_HPP
#define TINKER_DISTANCE_HPP

#include "utils/std.hpp"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"

class Frame;

class Distance : public AbstractAnalysis {
public:

    Distance();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Distance between two groups (mass-weighted)"; }

protected:

    template<typename SinglePassRange>
    [[nodiscard]] static std::tuple<double, double, double> calculate_mass_center(const SinglePassRange &atoms_group);

    std::deque<double> distances;

    Atom::AmberMask mask_for_group1;
    Atom::AmberMask mask_for_group2;

    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group1;
    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group2;

    boost::accumulators::accumulator_set<
            double,
            boost::accumulators::features<
                    boost::accumulators::tag::mean,
                    boost::accumulators::tag::variance
            >
    > acc;
};

#endif //TINKER_DISTANCE_HPP
