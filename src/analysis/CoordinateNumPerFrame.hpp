//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_COORDINATENUMPERFRAME_HPP
#define TINKER_COORDINATENUMPERFRAME_HPP

#include <unordered_set>
#include <memory>
#include <string>
#include <list>
#include "common.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class CoordinateNumPerFrame : public AbstractAnalysis {
public:

    CoordinateNumPerFrame();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Coordinate Number per Frame"; }

private:

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;


    double dist_cutoff;
    std::list<int> cn_list;
};


#endif //TINKER_COORDINATENUMPERFRAME_HPP
