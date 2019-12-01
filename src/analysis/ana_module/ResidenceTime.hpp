//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_RESIDENCETIME_HPP
#define TINKER_RESIDENCETIME_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <Eigen/Eigen>
#include <tbb/tbb.h>
#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"

class Frame;

class ResidenceTime : public AbstractAnalysis {
public:

    ResidenceTime();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void setParameters(const Atom::Node &id1, const Atom::Node &id2,
                       double cutoff, int t_star, const std::string &outfilename);

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "ResidenceTime"; }

protected:

    std::map<int, std::list<bool>> mark_map;

    void calculate();

    double dis_cutoff{};
    int steps = 0;
    int atom_num = 0;
    std::vector<double, tbb::tbb_allocator<double>> Rt_array;

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> mark;

    int time_star = 0;

    enum class TimeStarMode : int {
        Loose = 0,
        Strict = 1
    } timeStarMode;

    AmberMask center_atom_mask;
    AmberMask Ow_atom_mask;

    std::unordered_set<std::shared_ptr<Atom>> center_atom_group;
    std::unordered_set<std::shared_ptr<Atom>> Ow_atom_group;

    void fill_mark();

    void readTimeStarSetting();

    void readBasicSetting();
};

#endif //TINKER_RESIDENCETIME_HPP
