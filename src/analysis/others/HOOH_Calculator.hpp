

#ifndef TINKER_HOOH_CALCULATOR_HPP
#define TINKER_HOOH_CALCULATOR_HPP

#include <Eigen/Eigen>
#include <string_view>

class HOOH_Calculator {
   public:
    static void process();

    [[nodiscard]] static std::string_view title() { return "H2O2 coordinate calculator"; }

    static std::array<Eigen::Vector3d,2> calculate(Eigen::Vector3d &A, const Eigen::Vector3d &B, Eigen::Vector3d &C, double l_BA,
                                     double l_BC, double theta, double phi);
};

#endif  // TINKER_HOOH_CALCULATOR_HPP
