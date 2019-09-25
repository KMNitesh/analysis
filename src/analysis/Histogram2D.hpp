//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_HISTOGRAM2D_HPP
#define TINKER_HISTOGRAM2D_HPP

#include <utility>
#include <map>
#include <tuple>
#include <vector>

class Histogram2D {
public:
    Histogram2D(const std::pair<double, double> &range1, double width1,
                const std::pair<double, double> &range2, double width2);

    Histogram2D() = default;

    void update(double value1, double value2);
    void update(std::pair<double,double> value_pair);

    [[nodiscard]] std::vector<std::tuple<double, double, double>> getDistribution() const;

    [[nodiscard]] std::pair<double, double> getWidths() const { return {dimension1_width, dimension2_width}; }

    void initialize(const std::pair<double, double> &range1, double width1,
                    const std::pair<double, double> &range2, double width2);

    void initialize(double range1_max, double width1, double range2_max, double width2);

protected:

    std::pair<double, double> dimension1_range;
    std::pair<double, double> dimension2_range;
    double dimension1_width;
    double dimension2_width;
    int dimension1_bins;
    int dimension2_bins;

    std::map<std::pair<int, int>, size_t> hist;
};


#endif //TINKER_HISTOGRAM2D_HPP
