//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_VECTORSELECTOR_HPP
#define TINKER_VECTORSELECTOR_HPP

#include <tuple>
#include <vector>
#include <ostream>
#include <memory>

class Frame;

class Molecule;

class VectorSelector {
public:

    virtual int initialize(const std::shared_ptr<Frame> &frame) = 0;

    [[nodiscard]] virtual std::vector<std::tuple<double, double, double>>
    calculateVectors(const std::shared_ptr<Frame> &frame) = 0;

    [[nodiscard]] virtual std::tuple<double, double, double>
    calculateVector(const std::shared_ptr<Molecule> &mol, const std::shared_ptr<Frame> &frame) = 0;

    virtual void readInfo() = 0;

    [[nodiscard]] virtual std::string description() = 0;

    virtual void print(std::ostream &os) = 0;

    virtual ~VectorSelector() = default;
};

inline std::ostream &operator<<(std::ostream &os, std::shared_ptr<VectorSelector> &selector) {
    selector->print(os);
    return os;
}


#endif //TINKER_VECTORSELECTOR_HPP
