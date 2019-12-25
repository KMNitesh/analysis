#ifndef TINKER_GAUSSIANFCHKREADER_HPP
#define TINKER_GAUSSIANFCHKREADER_HPP

#include "utils/std.hpp"

class GaussianFchkReader {
public:

    static std::vector<double> readCartesian(std::istream &is);

    static std::tuple<double, double, double> getCoordinateOfAtom(std::size_t index /* low bound 1*/,
                                                                  const std::vector<double> &cartesians);

    void readCartesian(const std::string &filename) {
        std::ifstream ifstream(filename);
        cartesians = readCartesian(ifstream);
    }

    [[nodiscard]] std::tuple<double, double, double> getCoordinateOfAtom(std::size_t index /* low bound 1*/) const {
        return getCoordinateOfAtom(index, cartesians);
    }

    [[nodiscard]] std::size_t getTotalAtomNumbers() const { return cartesians.size() / 3; }

private:

    std::vector<double> cartesians;
};

#endif //TINKER_GAUSSIANFCHKREADER_HPP
