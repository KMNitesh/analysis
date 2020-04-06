//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_COORDINATENUMPERFRAME_HPP
#define TINKER_COORDINATENUMPERFRAME_HPP

#include <list>
#include <memory>
#include <string>
#include <unordered_set>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/common.hpp"

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
    AmberMask ids1;
    AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    double dist_cutoff;
    std::list<int> cn_list;
};

#endif  // TINKER_COORDINATENUMPERFRAME_HPP
