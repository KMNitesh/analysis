//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DISTANCE_HPP
#define TINKER_DISTANCE_HPP


#include <memory>
#include <list>
#include <string>
#include <unordered_set>
#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

// Distance
class Distance : public BasicAnalysis {
    std::list<double> group_dist_list;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;


public:
    Distance() {
        enable_outfile = true;
    }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const std::string title() {
        return "Distance";
    }

};

#endif //TINKER_DISTANCE_HPP
