//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_FINDMINBETWEENTWOGROUPS_HPP
#define TINKER_FINDMINBETWEENTWOGROUPS_HPP

#include <list>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <utility>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"

class Frame;
class Molecule;

class FindMinBetweenTwoGroups : public AbstractAnalysis {
public:
    FindMinBetweenTwoGroups();

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Find Min distance between two groups"; }

private:
    Atom::AmberMask amberMask;

    int total_frames = 0;

    std::vector<std::shared_ptr<Molecule>> mol_list;

    std::list<std::vector<double>> results;
};

#endif  // TINKER_FINDMINBETWEENTWOGROUPS_HPP
