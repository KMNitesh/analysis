//
// Created by xiamr on 7/4/19.
//

#include "Histogram.hpp"

using namespace std;

Histogram::Histogram(pair<double, double> dimension_range, double width) { initialize(dimension_range, width); }

void Histogram::initialize(double range_max, double width) { initialize({0, range_max}, width); }

void Histogram::initialize(const pair<double, double> &range, double width) {
    dimension_range = range;
    dimension_width = width;

    auto[dimension_min, dimension_max] = range;

    dimension_bins = static_cast<int>((dimension_max - dimension_min) / dimension_width);

    for (int ibin = 1; ibin <= dimension_bins; ibin++) {
        hist[ibin] = 0;
    }
}

void Histogram::update(double value) {
    int ibin = int((value - dimension_range.first) / dimension_width) + 1;

    if (ibin <= dimension_bins) {
        hist[ibin] += 1;
    }
}

std::vector<pair<double, double>> Histogram::getDistribution() const {
    double total = 0.0;

    for (int ibin = 1; ibin <= dimension_bins; ibin++) {
        total += hist.at(ibin);
    }

    std::vector<pair<double, double>> distribution;

    for (int ibin = 1; ibin <= dimension_bins; ibin++) {
        distribution.emplace_back((ibin - 0.5) * dimension_width + dimension_range.first,
                                  (hist.at(ibin) / total) / dimension_width);
    }
    return distribution;
}
