//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RESIDENCETIME_HPP
#define TINKER_RESIDENCETIME_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;

class ResidenceTime : public BasicAnalysis {

public:

    ResidenceTime() {
        enable_tbb = true;
        enable_outfile = true;
    }

    ~ResidenceTime();;

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void setParameters(const Atom::Node &id1, const Atom::Node &id2,
                       double cutoff, int t_star, const std::string &outfilename);

    void readInfo() override;

    static const std::string title() {
        return "ResidenceTime";
    }

private:

    std::map<int, std::list<bool>> mark_map;

    void calculate();

    void setSteps(size_t steps, int atom_num);

    double dis_cutoff;
    size_t steps = 0;
    int atom_num = 0;
    int *time_array = nullptr;
    double *Rt_array = nullptr;
    int **mark = nullptr;
    double time_star = 0;

    Atom::AtomIndenter ids1;
    Atom::AtomIndenter ids2;

    std::unordered_set<std::shared_ptr<Atom>> group1;
    std::unordered_set<std::shared_ptr<Atom>> group2;

};

#endif //TINKER_RESIDENCETIME_HPP
