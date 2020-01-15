#ifndef TINKER_GAUSSIANLOGCOORDINATEFORMAT_HPP
#define TINKER_GAUSSIANLOGCOORDINATEFORMAT_HPP

#include <string_view>

class GaussianLogCoordinateFormat {
public:
    [[nodiscard]] static std::string_view title() { return "Gaussian log Coordinate Format"; }

    static void process();
};


#endif //TINKER_GAUSSIANLOGCOORDINATEFORMAT_HPP
