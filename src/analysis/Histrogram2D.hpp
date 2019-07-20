//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_HISTROGRAM2D_HPP
#define TINKER_HISTROGRAM2D_HPP

#include <utility>
#include <map>
#include <tuple>
#include <vector>

class Histrogram2D {
public:
    Histrogram2D(const std::pair<double, double> &range1, double width1,
                 const std::pair<double, double> &range2, double width2);

    Histrogram2D() = default;

    void update(double value1, double value2);

    std::vector<std::tuple<double, double, double>> getDistribution() const;

    std::pair<double, double> getWidths() const { return {dimension1_width, dimension2_width}; }

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


#endif //TINKER_HISTROGRAM2D_HPP
