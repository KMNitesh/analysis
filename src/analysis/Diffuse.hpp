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
    Diffuse() {
        enable_forcefield = true;
        enable_outfile = true;
    }

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

    std::unordered_set<std::shared_ptr<Atom>> group;
    bool first_round = true;

    bool first_frame = true;
    int total_frame_number;
    int steps = 0;
    double time_increment_ps = 0.1;
    int total_mol = 0;


    bool bSerial = true;
    bool bTradition = true;

    Eigen::MatrixXd xcm;
    Eigen::MatrixXd ycm;
    Eigen::MatrixXd zcm;

};


#endif //TINKER_DIFFUSE_HPP
