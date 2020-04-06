//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIFFUSE_HPP
#define TINKER_DIFFUSE_HPP

#include <Eigen/Eigen>

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/common.hpp"
#include "utils/std.hpp"

class Frame;

class Diffuse : public AbstractAnalysis {
public:
    Diffuse();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] std::string description() override;

    void setParameters(const AmberMask &mask, double time_increment_ps, int total_frames,
                       const std::string &outfilename);

    [[nodiscard]] static std::string_view title() {
        return "Self-Diffusion Coefficient Calculation based on Einstein Equation";
    }

private:
    AmberMask mask;

    std::set<std::shared_ptr<Molecule>> mols;
    bool first_round = true;
    int total_frame_number;
    int steps = 0;
    double time_increment_ps = 0.1;

    Eigen::Matrix<std::tuple<double, double, double>, Eigen::Dynamic, Eigen::Dynamic> xyzcm;
};

#endif  // TINKER_DIFFUSE_HPP
