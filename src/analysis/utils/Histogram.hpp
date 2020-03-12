//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_HISTOGRAM_HPP
#define TINKER_HISTOGRAM_HPP

#include <map>
#include <utility>
#include <vector>

class Histogram {
public:
    Histogram(std::pair<double, double> dimension_range, double width);

    Histogram() = default;

    void update(double value);

    std::vector<std::pair<double, double>> getDistribution() const;

    double getWidth() const { return dimension_width; }

    void initialize(const std::pair<double, double> &range, double width);

    void initialize(double range_max, double width);

    // protected:
    std::pair<double, double> dimension_range;
    double dimension_width;
    int dimension_bins;

    std::map<int, size_t> hist;
};

#endif  // TINKER_HISTOGRAM_HPP
