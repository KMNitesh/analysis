//
// Created by xiamr on 7/24/19.
//

#ifndef TINKER_TINKERDYNREADER_HPP
#define TINKER_TINKERDYNREADER_HPP

#include <fstream>
#include <memory>
#include <tuple>
#include <vector>

struct PBCBoxStruct {
    double xbox, ybox, zbox;
    double alpha, beta, gamma;
};

class AtomicPosition {
public:
    std::vector<std::tuple<double, double, double>> coordinates;
};

class TinkerDynReader {
public:
    explicit TinkerDynReader(std::shared_ptr<std::istream> istream) : istream(std::move(istream)){};

    void readContent();

    const PBCBoxStruct &getPBCBox() const { return box; }

    const AtomicPosition &getAtomicPosition() const { return positions; }

protected:
    std::shared_ptr<std::istream> istream;
    PBCBoxStruct box{};
    AtomicPosition positions;
};

#endif  // TINKER_TINKERDYNREADER_HPP
