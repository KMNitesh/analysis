//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DISTANCE_HPP
#define TINKER_DISTANCE_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "utils/std.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class Distance : public AbstractAnalysis {
public:
    Distance();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Distance between two groups (mass-weighted)"; }

    void setParameters(const AmberMask &A, const AmberMask &B, const std::string &out);

protected:
    template <typename SinglePassRange>
    [[nodiscard]] static std::tuple<double, double, double> calculate_mass_center(const SinglePassRange &atoms_group);

    void saveJson(std::ostream &os) const;

    std::deque<double> distances;

    AmberMask mask_for_group1;
    AmberMask mask_for_group2;

    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group1;
    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group2;

    boost::accumulators::accumulator_set<
        double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>>
        acc;
};

#endif  // TINKER_DISTANCE_HPP
