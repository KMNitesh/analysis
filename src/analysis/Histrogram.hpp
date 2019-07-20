//
// Created by xiamr on 7/4/19.
//

#ifndef TINKER_HISTROGRAM_HPP
#define TINKER_HISTROGRAM_HPP

#include <utility>
#include <vector>
#include <map>

class Histrogram {
public:
    Histrogram(std::pair<double, double> dimension_range, double width);

    Histrogram() = default;

    void update(double value);

    std::vector<std::pair<double, double>> getDistribution() const;


    double getWidth() const { return dimension_width; }

    void initialize(const std::pair<double, double> &range, double width);

    void initialize(double range_max, double width);

//protected:
    std::pair<double, double> dimension_range;
    double dimension_width;
    int dimension_bins;

    std::map<int, size_t> hist;
};


#endif //TINKER_HISTROGRAM_HPP
