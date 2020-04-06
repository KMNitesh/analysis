
#ifndef CURVATURE_HPP
#define CURVATURE_HPP

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"
#include "utils/PBCUtils.hpp"

class Frame;

class Curvature : public AbstractAnalysis {
public:
    Curvature();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Curvature"; }

private:
    AmberMask mask;
    std::vector<std::shared_ptr<Atom>> atoms;

    PBCUtils::MolPair mol;

    double calculate_curvature();

    std::deque<double> curvature_lists;
};

#endif // CURVATURE_HPP