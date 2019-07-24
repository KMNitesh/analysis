//
// Created by xiamr on 7/24/19.
//

#ifndef TINKER_TINKERDYNREADER_HPP
#define TINKER_TINKERDYNREADER_HPP

#include <memory>
#include <fstream>
#include <tuple>
#include <vector>

class PBCBox {
public:
    double xbox, ybox, zbox;
    double alpha, beta, gamma;
};

class AtomicPosition {
public:
    std::vector<std::tuple<double, double, double>> coordinates;
};

class TinkerDynReader {
public:
    explicit TinkerDynReader(std::shared_ptr<std::istream> istream) : istream(std::move(istream)) {};

    void readContent();

    const PBCBox &getPBCBox() const { return box; }

    const AtomicPosition &getAtomicPosition() const { return positions; }

protected:
    std::shared_ptr<std::istream> istream;
    PBCBox box{};
    AtomicPosition positions;
};


#endif //TINKER_TINKERDYNREADER_HPP
