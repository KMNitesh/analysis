//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DEMIXINDEXOFTWOGROUP_HPP
#define TINKER_DEMIXINDEXOFTWOGROUP_HPP


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

class DemixIndexOfTwoGroup : public BasicAnalysis {
public:
    DemixIndexOfTwoGroup() { enable_outfile = true; }

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void setParameters(const Atom::Node &id1, const Atom::Node &id2, const Grid &grid, const std::string &outfilename);

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    std::string getOutfileName() override { return outfilename; }

    static const std::string title() { return "Calculate demix index of two groups"; }

private:

    auto calculate_grid_index(const std::shared_ptr<Atom> &atom, const std::shared_ptr<Frame> &frame);

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    int grid_x;
    int grid_y;
    int grid_z;

    std::list<std::tuple<double, double>> demix_index_list;

    std::string outfilename;
};


#endif //TINKER_DEMIXINDEXOFTWOGROUP_HPP
