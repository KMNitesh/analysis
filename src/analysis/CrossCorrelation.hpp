//
// Created by xiamr on 8/28/19.
//

#ifndef TINKER_CROSSCORRELATION_HPP
#define TINKER_CROSSCORRELATION_HPP

#include "std.hpp"

class CrossCorrelation {
public:
    static std::string title() { return "Cross-correlation Analysis"; }

    static void calculate(const std::string &out);

    static double calCovariance(const std::deque<std::tuple<double, double, double>> &series);

    static std::vector<double> calculateCrossCorrelation(std::deque<std::tuple<double, double, double>> &series);

    static std::vector<double> calculateConvolutionFunction(const std::deque<std::tuple<double, double, double>> &series);
};


#endif //TINKER_CROSSCORRELATION_HPP
