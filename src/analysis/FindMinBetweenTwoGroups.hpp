//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_FINDMINBETWEENTWOGROUPS_HPP
#define TINKER_FINDMINBETWEENTWOGROUPS_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <utility>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class Molecule;


class FindMinBetweenTwoGroups : public BasicAnalysis {
public:
    FindMinBetweenTwoGroups() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    static const std::string title() { return "Find Min distance between two groups"; }

private:
    Atom::AmberMask ids;

    int total_frames = 0;

    std::vector<std::shared_ptr<Molecule>> mol_list;

    std::list<std::vector<double>> results;
};

#endif //TINKER_FINDMINBETWEENTWOGROUPS_HPP
