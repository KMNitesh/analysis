//
// Created by xiamr on 8/6/19.
//

#ifndef TINKER_HBONDLIFETIME_HPP
#define TINKER_HBONDLIFETIME_HPP

#include "std.hpp"
#include "BasicAnalysis.hpp"

class Frame;

class Atom;

class HBondLifeTime : public BasicAnalysis {
public:
    HBondLifeTime();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() { return "Hydrogen Bond LifeTime for Water(Ci, history independent)"; }

protected:

    double dist_R_cutoff;
    double angle_HOO_cutoff;

    double time_increment_ps;
    double max_time_grap_ps;

    std::vector<std::deque<int>> hb_histroy;

    std::deque<std::shared_ptr<Atom>> hydrogens;
    std::deque<std::shared_ptr<Atom>> oxygens;

    std::vector<double> calculateAcf() const;

    void printData(std::ostream &os, const std::vector<double> &acf) const;
};


#endif //TINKER_HBONDLIFETIME_HPP
