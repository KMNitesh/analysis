//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIFFUSE_HPP
#define TINKER_DIFFUSE_HPP

#include "std.hpp"

#include <Eigen/Eigen>

#include "common.hpp"
#include "AbstractAnalysis.hpp"
#include "atom.hpp"

class Frame;


class Diffuse : public AbstractAnalysis {

public:
    Diffuse();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    std::string description() override;

    void
    setParameters(const Atom::Node &mask, double time_increment_ps, int total_frames, const std::string &outfilename);

    static const std::string title() {
        return "Self-Diffusion Coefficient Calculation based on Einstein Equation";
    }

private:

    Atom::AmberMask ids;

    std::set<std::shared_ptr<Molecule>> mols;
    bool first_round = true;
    int total_frame_number;
    int steps = 0;
    double time_increment_ps = 0.1;

    Eigen::Matrix<Eigen::Array3d, Eigen::Dynamic, Eigen::Dynamic> xyzcm;

};


#endif //TINKER_DIFFUSE_HPP
