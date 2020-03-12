//
// Created by xiamr on 8/28/19.
//

#ifndef TINKER_CROSSCORRELATION_HPP
#define TINKER_CROSSCORRELATION_HPP

#include "utils/std.hpp"

class CrossCorrelation {
public:
    [[nodiscard]] static std::string_view title() { return "Cross-correlation Analysis"; }

    static void calculate(const std::string &out);

    [[nodiscard]] static double calCovariance(const std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static std::vector<double> calculateCrossCorrelation(
        std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static std::vector<double> calculateCrossCorrelation2(
        std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static std::vector<double> calculateConvolutionFunction(
        const std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static std::deque<std::tuple<double, double, double>> readAngleSeries();

    [[nodiscard]] static std::vector<std::tuple<double, double, double>> getDistribution(
        const std::deque<std::tuple<double, double, double>> &series);

    static void printCrossCorrelationFunction(const std::vector<double> &cross_correlation_function,
                                              std::deque<std::tuple<double, double, double>> &series,
                                              std::ofstream &os);

    static void printCrossCorrelationFunction2(const std::vector<double> &cross_correlation_function,
                                               std::deque<std::tuple<double, double, double>> &series,
                                               std::ofstream &os);

    static void printConvolutionFunction(const std::vector<double> &convolution_function,
                                         std::deque<std::tuple<double, double, double>> &series, std::ofstream &os);

    static void printHistrogramDistribution(const std::vector<std::tuple<double, double, double>> &distribution,
                                            std::ofstream &os);

    [[nodiscard]] static std::pair<double, double> calcSeriesAverage(
        const std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static double calculateCrossCorrelationDWang(std::deque<std::tuple<double, double, double>> &series);

    [[nodiscard]] static double calculateCrossCorrelationDWang2(std::deque<std::tuple<double, double, double>> &series);
};

#endif  // TINKER_CROSSCORRELATION_HPP
