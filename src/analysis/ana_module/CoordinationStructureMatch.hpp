//
// Created by xiamr on 10/12/19.
//

#ifndef TINKER_COORDINATIONSTRUCTUREMATCH_HPP
#define TINKER_COORDINATIONSTRUCTUREMATCH_HPP

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "utils/std.hpp"

class Frame;

class CoordinationStructureMatch : public AbstractAnalysis {
public:
    CoordinationStructureMatch();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() {
        return "Coordination Structure Match based on Nonlinear Least Squares Pattern Recognition";
    }

    [[nodiscard]] static double testCASP(std::vector<std::tuple<double, double, double>> &coord);

    [[nodiscard]] static double testTCTP(std::vector<std::tuple<double, double, double>> &coord);

protected:
    AmberMask metal_mask;
    AmberMask Ow_atom_mask;

    std::shared_ptr<Atom> metal;
    std::vector<std::shared_ptr<Atom>> Ow_atoms;
    double cutoff2;

    // pair<TCTP,CASP>
    std::deque<std::pair<double, double>> r_list;
};

#endif  // TINKER_COORDINATIONSTRUCTUREMATCH_HPP
