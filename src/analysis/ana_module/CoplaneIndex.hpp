
#ifndef TINKER_COPLANEINDEX_HPP
#define TINKER_COPLANEINDEX_HPP

#include "AbstractAnalysis.hpp"
#include "data_structure/atom.hpp"
#include "dsl/AmberMask.hpp"

class Frame;

class CoplaneIndex : public AbstractAnalysis {
public:
    CoplaneIndex();

    void processFirstFrame(std::shared_ptr<Frame> &frame) override;

    void process(std::shared_ptr<Frame> &frame) override;

    void print(std::ostream &os) override;

    void readInfo() override;

    [[nodiscard]] static std::string_view title() { return "Coplane index"; }

private:
    std::vector<std::array<AmberMask, 3>> mask_arrays;

    std::vector<std::array<std::shared_ptr<Atom>, 3>> atom_arrays;

    std::deque<std::pair<double, double>> coplaneIndex;  // mean and std
};

#endif  // TINKER_COPLANEINDEX_HPP
