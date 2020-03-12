//
// Created by xiamr on 7/4/19.
//

#include "Histogram2D.hpp"

Histogram2D::Histogram2D(const std::pair<double, double> &range1, double width1,
                         const std::pair<double, double> &range2, double width2) {
    initialize(range1, width1, range2, width2);
}

void Histogram2D::update(double value1, double value2) {
    int ibin1 = int((value1 - dimension1_range.first) / dimension1_width) + 1;
    int ibin2 = int((value2 - dimension2_range.first) / dimension2_width) + 1;

    if (ibin1 <= dimension1_bins and ibin2 <= dimension2_bins) {
        hist[{ibin1, ibin2}] += 1;
    }
}

void Histogram2D::update(std::pair<double, double> value_pair) {
    auto [value1, value2] = value_pair;
    update(value1, value2);
}

std::vector<std::tuple<double, double, double>> Histogram2D::getDistribution() const {
    std::vector<std::tuple<double, double, double>> distribution;
    for (int ibin1 = 1; ibin1 <= dimension1_bins; ibin1++) {
        for (int ibin2 = 1; ibin2 <= dimension2_bins; ibin2++) {
            distribution.emplace_back((ibin1 - 0.5) * dimension1_width + dimension1_range.first,
                                      (ibin2 - 0.5) * dimension2_width + dimension2_range.first,
                                      hist.at({ibin1, ibin2}) / (dimension1_width * dimension2_width));
        }
    }
    return distribution;
}

void Histogram2D::initialize(const std::pair<double, double> &range1, double width1,
                             const std::pair<double, double> &range2, double width2) {
    dimension1_range = range1;
    dimension1_width = width1;

    dimension2_range = range2;
    dimension2_width = width2;

    auto [dimension1_min, dimension1_max] = range1;
    auto [dimension2_min, dimension2_max] = range2;

    dimension1_bins = static_cast<int>((dimension1_max - dimension1_min) / dimension1_width);
    dimension2_bins = static_cast<int>((dimension2_max - dimension2_min) / dimension2_width);

    for (int ibin1 = 1; ibin1 <= dimension1_bins; ibin1++) {
        for (int ibin2 = 1; ibin2 <= dimension2_bins; ibin2++) {
            hist[{ibin1, ibin2}] = 0;
        }
    }
}

void Histogram2D::initialize(double range1_max, double width1, double range2_max, double width2) {
    initialize({0.0, range1_max}, width1, {0.0, range2_max}, width2);
}
