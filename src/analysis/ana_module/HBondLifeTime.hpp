//
// Created by xiamr on 8/6/19.
//

#ifndef TINKER_HBONDLIFETIME_HPP
#define TINKER_HBONDLIFETIME_HPP

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/std.hpp"

class Frame;

class HBondLifeTime : public AbstractAnalysis {
public:
    HBondLifeTime();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() {
        return "Hydrogen Bond LifeTime for Water(Ci, history independent)";
    }

protected:
    AmberMask water_mask;

    double dist_R_cutoff;
    double angle_HOO_cutoff;

    double time_increment_ps;
    double max_time_grap_ps;

    std::vector<std::deque<int>> hb_histroy;

    std::deque<std::shared_ptr<Atom>> hydrogens;
    std::deque<std::shared_ptr<Atom>> oxygens;

    [[nodiscard]] virtual std::vector<double> calculateAcf() const;

    void printData(std::ostream &os, const std::vector<double> &acf, std::string_view title) const;
};

#endif  // TINKER_HBONDLIFETIME_HPP
