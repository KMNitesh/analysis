//
// Created by xiamr on 6/14/19.
//

#ifndef TINKER_DIFFUSE_HPP
#define TINKER_DIFFUSE_HPP

#include <memory>
#include <unordered_set>
#include <string>
#include <map>
#include <list>
#include <utility>

#include <Eigen/Eigen>

#include "common.hpp"
#include "BasicAnalysis.hpp"
#include "atom.hpp"

class Frame;


class Diffuse : public BasicAnalysis {

public:
    Diffuse() {
        enable_forcefield = true;
        enable_outfile = true;
    }

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    static const std::string title() {
        return "Diffusion Coefficient by Einstein equation";
    }

private:

    Atom::AtomIndenter ids;

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
