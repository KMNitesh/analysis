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
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;

class DemixIndexOfTwoGroup : public AbstractAnalysis {
public:

    DemixIndexOfTwoGroup();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    void setParameters(const Atom::Node &id1, const Atom::Node &id2, const Grid &grid, const std::string &outfilename);


    [[nodiscard]] static std::string_view title() { return "Calculate demix index of two groups"; }

private:

    [[nodiscard]] auto calculate_grid_index(const std::shared_ptr<Atom> &atom, const std::shared_ptr<Frame> &frame);

    Atom::AmberMask ids1;
    Atom::AmberMask ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

    int grid_x;
    int grid_y;
    int grid_z;

    std::list<std::tuple<double, double>> demix_index_list;

};


#endif //TINKER_DEMIXINDEXOFTWOGROUP_HPP
