//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DISTANCE_HPP
#define TINKER_DISTANCE_HPP

#include "std.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

// Distance
class Distance : public AbstractAnalysis {

public:
    Distance();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static std::string title() { return "Distance between two groups (mass-weighted)"; }

protected:

    template<typename SinglePassRange>
    std::tuple<double, double, double> calculate_mass_center(SinglePassRange &&atoms_group);

    std::deque<double> distances;

    Atom::AmberMask mask_for_group1;
    Atom::AmberMask mask_for_group2;

    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group1;
    std::unordered_set<std::shared_ptr<Atom>> atoms_for_group2;
};

#endif //TINKER_DISTANCE_HPP
