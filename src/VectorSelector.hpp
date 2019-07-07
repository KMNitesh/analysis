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

    virtual std::vector<std::tuple<double, double, double>> calculateVectors(const std::shared_ptr<Frame> &frame) = 0;

    virtual std::tuple<double, double, double> calculateVector(const std::shared_ptr<Molecule> &mol,
                                                               const std::shared_ptr<Frame> &frame) = 0;

    virtual void readInfo() = 0;

    virtual void print(std::ostream &os) = 0;

    virtual ~VectorSelector() = default;
};


#endif //TINKER_VECTORSELECTOR_HPP
